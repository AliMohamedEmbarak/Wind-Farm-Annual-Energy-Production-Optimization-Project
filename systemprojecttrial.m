clc
clear variables

%% Turbine Specifications.
turb_specs        = struct();
turb_specs.name   = "Anon";
turb_specs.vendor = "Anon";
turb_specs.type   = "Anon";
turb_specs.diameter  = 100.0;  % in m
turb_specs.rotorArea = 7853; % in m2
turb_specs.hubHeight = 100;  % in m
turb_specs.cutInWindSpeed  = 3.5; % in m/s
turb_specs.cutOutWindSpeed = 25;  % in m/s
turb_specs.ratedWindSpeed  = 15;  % in m/s
turb_specs.ratedPower = 3; % in MW

wind_drct= 90;
wind_sped= 30;

nTurbines = 25;
farmSize = 4000; % in meters
clearance = 50; % in meters
minDistance = 400; % in meters


%% Checking if the turbine coordinates satisfies the constraints
% Comment out the function call to checkConstraints below if you desire. 
% Note that this is just a check and the function does not quantifies the 
% amount by which the constraints are violated if any. 
%checkConstraints(positions, turb_specs.diameter)

%% Load Data Files and Read
% Turbine x,y coordinates
layoutFile = 'Shell_Hackathon Dataset\turbine_loc_test.csv';
% Power curve
powerCurveFile = 'Shell_Hackathon Dataset\power_curve.csv';
% Wind Data
windDataFile = 'Shell_Hackathon Dataset\wind_data\wind_data_2007.csv';

%turb_coords  = table2array(readtable(layoutFile));    
power_curve  = readtable(powerCurveFile);             
windData     = readtable(windDataFile);                

%% 
% Bin wind data and calculate probability of wind instance occurence
% -**-SHOULD NOT BE MODIFIED-**-
% speed 'slices'
% n_slices_sped = 15;
% slices_sped   = linspace(0,30.0,n_slices_sped+1);
% 
% % direction 'slices' in degrees
% slices_drct = [360.0,linspace(10.0,350.0,35)];
% n_slices_drct     = length(slices_drct);
% 
% binned_wind = zeros(n_slices_drct,n_slices_sped);
% 
% for i = 1:n_slices_drct
%     for j = 1:n_slices_sped
%         foo = windData.drct == slices_drct(i) & ...
%                windData.sped >= slices_sped(j) &  ...
%                windData.sped < slices_sped(j+1);
%         
%         binned_wind(i,j) = sum(foo);         
%     end
% end
% wind_inst_freq = binned_wind/sum(sum(binned_wind));
% directionBins  = slices_drct;
% % take the mid value as effective speed
% speedBins      = mean([slices_sped(1:end-1);slices_sped(2:end)]);
  

%% calling functions
positions = initialvalues(nTurbines, farmSize, clearance, minDistance);
power = partAEP(positions,power_curve,wind_drct,wind_sped,turb_specs);

%% AEP calculation
disp('Calculating AEP......');
[AEP] = totalAEP(power);
format long
disp(AEP)
%%
figure
plot(positions(:,1),positions(:,2),'ro');
title(['AEP = ', num2str(AEP),' GW.hr'])
xlabel('Farm Width [m]');
ylabel('Farm length [m]')

%% clear temporary vars
% clearvars -except AEP windData windDataBinCount WindTurbineParameters 

%% Initial positions function
function [positions] = initialvalues(nTurbines, farmSize, clearance, minDistance)
    
    % Calculate the available space for turbine placement
    availableSize = farmSize - 2 * clearance;

    % Check if the available space is sufficient
    if availableSize <= 0
        error("Insufficient space for turbine placement");
    end

    % Calculate the grid size for turbine placement
    gridSize = availableSize / (sqrt(nTurbines) - 1);

    % Initialize position variable
    positions = zeros(2 * nTurbines, 1);

    % Adjust positions to satisfy minimum distance constraint
    for i = 1:nTurbines
        x = clearance + mod(i-1, sqrt(nTurbines)) * gridSize;
        y = clearance + floor((i-1) / sqrt(nTurbines)) * gridSize;

        for j = 1:i-1
            dist = sqrt((x-positions((j-1)*2+1))^2 + (y-positions(j*2))^2);
            while dist < minDistance
                % Adjust the position of turbine i
                angle = rand * 2 * pi; % Random angle
                x = positions((j-1)*2+1) + minDistance * cos(angle);
                y = positions(j*2) + minDistance * sin(angle);
                dist = sqrt((x-positions((j-1)*2+1))^2 + (y-positions(j*2))^2);
            end
        end

        % Assign the position of turbine i
        positions((i-1)*2+1) = x;
        positions(i*2) = y;
        
        positions = transpose(positions);
    end
end

%%
function C = constraint(nTurbines)
  
    C = sym(zeros(1, 300));
    index = 1;
    
    x = sym('x', [1 nTurbines]);
    y = sym('y', [1 nTurbines]);
    
    for i = 1:nTurbines-1
        for j = i+1:nTurbines
            dist = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);
            C(index) = 400 - dist;
            index = index + 1;
        end
    end
    C = transpose(C);  % Transpose the constraints to one column
disp(C)
end

%% AEP calculation
%% 
% -**-SHOULD NOT BE MODIFIED-**-   
% Returns total power produced by all the turbine for this particular 
% 'wind instance' represented by wind_drct, wind_sped.

function power = partAEP(positions,power_curve,wind_drct,wind_sped,turb_specs)
% number of turbines in the layout
n_turbs = length(positions);

% Turbine radius
turb_rad = turb_specs.diameter/2.0;

% we use power_curve data as look up to estimate the thrust coeff.
% of the turbine for the corresponding closest matching wind speed
% Ct = interpn(power_curve.WindSpeed_m_s_,power_curve.ThrustCoeffecient,wind_sped,'nearest');
Ct = qinterp1(power_curve.WindSpeed_m_s_,power_curve.ThrustCoeffecient,wind_sped,0);

% Wake decay constant kw for the offshore case
kw = 0.05;

% converting wind direction from degrees to radians
wind_drct = wind_drct - 90.0;
wind_drct = deg2rad(wind_drct);

% Constants to be used in converting turbine eucledian coordinates to downwind
% and crosswind coordinates
cos_dir = cos(wind_drct);
sin_dir = sin(wind_drct);

% initializing turbine power vector
turb_pwr = zeros(n_turbs,1);

% Calculating wake effect of turbines on each other
impact_on_ibyj = zeros(n_turbs,n_turbs);

% i - target turbine
for i=1:n_turbs
    % looping over all other turbs to check their effect
    for j=1:n_turbs
        
        % Turbines seperation distance in downwind direction
        x = (positions(i,1) - positions(j,1))* cos_dir -  ...
            (positions(i,2) - positions(j,2)) * sin_dir;
        
        % Turbines seperation distance in crosswind direction
        y = (positions(i,1) - positions(j,1))* sin_dir + ...
            (positions(i,2) - positions(j,2)) * cos_dir;
        
        % Naturally, no wake effect of turbine on itself
        if i ~= j
            
            % either j not an upstream turbine or wake not happening 
            % on i because its outside of the wake region of j
            if (x <= 0) || (abs(y) > (turb_rad + kw*x))
                impact_on_ibyj(i,j) = 0.0;
            
            % otherwise, at target i, wake is happening due to j
            else
                impact_on_ibyj(i,j) = (1-sqrt(1-Ct)) * ...
                    ((turb_rad/(turb_rad + kw*x))^2);
            end
        end
    end
    
    % Calculate Total vel deficit from all upstream turbs, using sqrt of sum of sqrs
    sped_deficit = sqrt(sum(impact_on_ibyj(i,:) .* impact_on_ibyj(i,:)));
    
    % Effective wind speed at ith turbibe
    wind_sped_eff = wind_sped * (1.0 - sped_deficit);
    
    % we use power_curve data as look up to estimate the power produced
    % by the turbine for the corresponding closest matching wind speed
    % turb_pwr(i) = interpn(power_curve.WindSpeed_m_s_, power_curve.Power_MW_, ...
    %                           wind_sped_eff, 'nearest');
    turb_pwr(i) = qinterp1(power_curve.WindSpeed_m_s_, power_curve.Power_MW_, ...
                              wind_sped_eff, 0);                      
end
% Sum the power from all turbines for this wind instance
power = sum(turb_pwr);
end

function [AEP] = totalAEP(power)

% AEP in GWh
hrs_per_year = 365.0 * 24.0;
AEP = hrs_per_year * power/1000;

end

%%
function checkConstraints(turb_coords, turb_dia)

% -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-

% Checks if the turbine configuration satisfies the two
% constraints:(i) perimeter constraint,(ii) proximity constraint 
% Prints which constraints are violated if any. Note that this 
% function does not quantifies the amount by which the constraints 
% are violated if any. 

prox_constr_viol = 0;
peri_constr_viol = 0;

% Defining wind farm layout geometry
gridSize = [4000 4000]; % Rectangular grid
numberOfTurbines = 25; % Number of turbines
minAllowableDistance =  4*turb_dia; % in m
boundaryMargin = 50; % in m

xmin = 0+boundaryMargin;
xmax = gridSize(1) - boundaryMargin;

ymin = 0+boundaryMargin;
ymax = gridSize(2) - boundaryMargin;

xboundary = [xmin,xmin,xmax,xmax];
yboundary = [ymin,ymax,ymax,ymin];

% checks if for every turbine perimeter constraint is satisfied.
% breaks out if False anywhere
for i = 1:numberOfTurbines
    xq = turb_coords(i,1);
    yq = turb_coords(i,2);
    incheck = inpolygon(xq,yq,xboundary,yboundary);
    if (incheck == 0)
        peri_constr_viol = 1;
        break
    end
end

% checks if for every turbines proximity constraint is satisfied. 
% breaks out if False anywhere
for i = 1:numberOfTurbines
    for j = 1:numberOfTurbines
        
        if (i ~= j)
            xi = positions(i,1);
            yi = positions(i,2);
            xj = positions(j,1);
            yj = positions(j,2);
            
            dist_ij = sqrt((xi-xj)^2 + (yi-yj)^2);
            if (dist_ij < minAllowableDistance)
                prox_constr_viol = 1;
                break
            end
            
        end
    end    
end 

% print messages
if  (peri_constr_viol  == 1  && prox_constr_viol == 1)
    fprintf('Somewhere both perimeter constraint and proximity constraint are violated\n')
elseif (peri_constr_viol == 1  && prox_constr_viol == 0)
    fprintf('Somewhere perimeter constraint is violated\n')
elseif (peri_constr_viol == 0 && prox_constr_viol == 1)
    fprintf('Somewhere proximity constraint is violated\n')
else
    fprintf('Both perimeter and proximity constraints are satisfied !!\n')
end

end




%%
function Yi = qinterp1(x,Y,xi,methodflag)
% -**-THIS FUNCTION SHOULD NOT BE MODIFIED-**-
% Performs fast interpolation compared to interp1
% Original Source: https://tinyurl.com/y2adgpzk

% Forces vectors to be columns
x = x(:); xi = xi(:);
sx = size(x); sY = size(Y);
if sx(1)~=sY(1)
    if sx(1)==sY(2)
        Y = Y';
    else
        error('x and Y must have the same number of rows');
    end
end

if nargin>=4
    method=methodflag;
else
    method = 1;    % choose nearest-lower-neighbor, linear, etc.
                   % uses integer over string for speed
end

% Gets the x spacing
ndx = 1/(x(2)-x(1)); % one over to perform divide only once
xi = xi - x(1);      % subtract minimum of x

% Fills Yi with NaNs
s = size(Y);
if length(s)>2
    error('Y may only be one- or two-dimensional');
end
Yi = NaN*ones(length(xi),s(2));

switch method
    case 0 %nearest-neighbor method
        rxi = round(xi*ndx)+1;        % indices of nearest-neighbors
        flag = rxi<1 | rxi>length(x) | isnan(xi);
                                      % finds indices out of bounds
        nflag = ~flag;                % finds indices in bounds
        Yi(nflag,:) = Y(rxi(nflag),:);
    case 1 %linear interpolation method
        fxi = floor(xi*ndx)+1;          % indices of nearest-lower-neighbors
        flag = fxi<1 | fxi>length(x)-1 | isnan(xi);
                                        % finds indices out of bounds
        nflag = ~flag;                  % finds indices in bounds
        Yi(nflag,:) = (fxi(nflag)-xi(nflag)*ndx).*Y(fxi(nflag),:)+...
            (1-fxi(nflag)+xi(nflag)*ndx).*Y(fxi(nflag)+1,:);
end
end


%%fmincon

% options = optimoptions('fmincon','Display','iter');
% [x,fval] = fmincon(fun,X0,[],[],[],[],lb,ub,nonlcon,options);

%fun: AEP, X0: initialvalues (positions is a column vector of 50 elements (first 25= xcoordinates),(second 25 = ycoordinates),  
% lb =lb= ones(1, 50)*50;  ub = ones(1,50)*50; gives them in one row
% lb_transpose = lb.'; ub_transpose = ub.';  
%nonlcnonlcon: constraint given in one row
% C = transpose(C);  % Transpose the constraints to one column




%% ga function 


% options = optimoptions('ga','PlotFcn', {@gaplotbestf, @gaplotbestindiv},'UseParallel',true,'FunctionTolerance',1e-10,'InitialPopulationMatrix',X0);
% [x,fval,exitflag,output] = ga(fun,nvars,[],[],[],[],lb,ub,nonlcon,options);

%nvars = 50 


% % Patternsearch
% options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf,'UseParallel',
% true);
% [x,fval] = patternsearch(fun,X0,[],[],[],[],lb,ub,nonlcon,options);

