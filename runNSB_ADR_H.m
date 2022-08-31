%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test of two platelet species in the H domain
% NSB:  
%       rho[ u_t + u u_x + v u_y] =  - p_x + mu ( u_xx + u_yy) - mu alpha
%       rho[ v_t + u v_x + v v_y] =  - p_y + mu ( v_xx + v_yy) - mu alpha
%                          div(u) = 0
% Platelets:
%       (P^m)_t = - div[ W (uP^m - D grad(P^m) ] - k_adh (Pmax - P^b)P^m
%       (P^b)_t = k_adh (Pmax - P^b)P^m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;

% plotFlag:  1 -> FV Grid, 
%            2 -> pcolor during time stepping
%            3 -> Plots Stream Lines
% gridFlag:  0 -> uniform grid
%            1 -> nonuniform grid
% incr:      # of iterations between plots
    
gridFlag = 1;
plotFlag = 3;
incr = 50;
saveFlag = 1;
tol = 1e-16; % set the breaking point for max|u^{n+1} - u^n| < tol
Pressure_fix = 0; % Best to not fix pressure 

% Domain H within [0, Wleft +Wmid +Wright] x [0, Hlow +Hmid +Hup]
% All lengths are in cm
Wleft = 1e-2; Wmid = 1.5e-2; Wright = 1e-2; 
Hlow  = .75e-2; Hmid = .2e-2; Hup   = .75e-2;

% Set Nx, Ny and dt for given test --------------------------------------
Nx = 200; % For convergence test Nx = 32, 64,..., 512
Ny = 200;

% Time
Tf = 10; % s, maximum final time in seconds

% Physical Parameters Fluid 
mu = 2.62507e-2; % g / cm /s
rho = 1; % g / cm^3

% Volumetric Flow Rates at Upper Left and Upper Right Boundaries
Q1 = 2.77778e-2; % cm^2 / s
Q3 = 1.55556e-2; % cm^2 / s
QM = 4.25e-3;    % cm^2 / s

% Initial inlet/outlet pressures (from HCA) 10 * g/cm/s^2 = 1 Pa 
Pw_in = 10.3471e1;  % g/cm/s^2 : Pressure Inlet (wash)
Pw_out = 0;         % g/cm/s^2 : Pressure Outlet (wash)
Pb_in = 294.076e1;  % g/cm/s^2 :Pressure Inlet (blood)
Pb_out = 284.213e1; % g/cm/s^2 : Pressure Outlet (blood)

% % Volumetric Flow Rates at Upper Left and Upper Right Boundaries
% Q3 = 2.77778e-2; % cm^2 / s
% Q1 = 1.55556e-2; % cm^2 / s
% QM = 4.25e-3;    % cm^2 / s
% 
% % Initial inlet/outlet pressures (from HCA) 10 * g/cm/s^2 = 1 Pa 
% Pb_in = 10.3471e1;  % g/cm/s^2 : Pressure Inlet (wash)
% Pb_out = 0;         % g/cm/s^2 : Pressure Outlet (wash)
% Pw_in = 294.076e1;  % g/cm/s^2 :Pressure Inlet (blood)
% Pw_out = 284.213e1; % g/cm/s^2 : Pressure Outlet (blood)


% Characteristic Scales (fluid)
Xchar = 3e-4; % cm Size of a platelet
Uchar = 1.55556; % cm/s
Tchar = Xchar/Uchar; % s
Prchar = rho * Uchar^2;
Bchar = rho * Uchar / (mu * Xchar);
% Bchar = rho / mu / Xchar;
Re = rho * Uchar * Xchar / mu;

% Physical Params (Platelets)
Dp = 2.5e-7; % cm^2 / s 
Pmax = 6.67e10; % platelets per cm^3

% Characteristic Scales (Platelet Densities ADR)
Plchar = 2.5e8; % platelets/cm^3 or 250,000 / mm^3
Rchar = Uchar / Plchar / Xchar; % 1/cm^3/s 
Pe = Xchar * Uchar / Dp;


% Pack various input parmeters for passing to NS solver
flags = {plotFlag, gridFlag, incr, saveFlag, tol}; % Store flags in cell array
prms = [mu, rho, Re, Uchar, Xchar, Tchar, Prchar, Bchar, ...
        Plchar, Dp, Pe, Rchar,Pmax]; % Store Parameters in vector
P_io = [Pw_in, Pw_out, Pb_in, Pb_out, Pressure_fix]; % Store Parameters in vector
Q = [Q1, Q3];

D = [Nx,Ny,Wleft,Wmid,Wright,Hlow,Hmid,Hup]; % Domain Params.

[S,B, grd, N, Tf] = SolveNSB_ADR(D,Tf,prms,Q,P_io,flags);


dataName = ['data/NSB_ADRdata',num2str(gridFlag),'Nx',num2str(Nx),...
    'Tf',num2str(Tf),'Pfix',num2str(Pressure_fix),'.mat'];
save(dataName)


