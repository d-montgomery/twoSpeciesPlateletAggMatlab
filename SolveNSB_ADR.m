function [S,B, grd, N, Tf] = SolveNSB_ADR(D,Tfdim,prms,Q,P_io,flags)
% Solve rho[ u_t + u u_x + v u_y] = -p_x + mu ( u_xx + u_yy) + fx
%       rho[ v_t + u v_x + v v_y] = -p_y + mu ( v_xx + v_yy) + fy
% in the H domain given by parameters Hlow, Hmid, Hup, WL, Wmid, WR.
% Inputs:
% D    : Domain parameters D = [Nx,Ny,Wleft,Wmid,Wright,Hlow,Hmid,Hup]
% t    : is a row vector for time 0 < t <= Tf
% prms : a vector with parameters prms = [mu, rho, w, beta, h]
% IC   : a 3x1 cell array where IC{1} = ICu = @(x,y) ...;  for u - eqt
%                               IC{2} = ICv = @(x,y) ...;  for v - eqt.
%                               IC{3} = ICp = @(x,y) ...;  for pressure.
% u_coeffs : matrix containing BC coefficients for u-momentum eqt, 
%            associated with \partial \Omega_i, u_coeffs = [a5 a6; b5 b6];   
% v_coeffs : matrix containing BC coefficients for v-momentum eqt, 
%                                   v_coeffs = [c1 c2 c3 c4; d1 d2 d3 d4];
% BCu, BCv : cell arrays containing RHS functions for each BC.
%                                   BCu = {g1,  g2,  g3,  g4,  g5,  g6;
%                                         wu1, wu2, wu3, wu4, wu5, wu6};      
%                                   BCv = {h1,  h2,  h3,  h4,  h5,  h6;
%                                         wv1, wv2, wv3, wv4, wv5, wv6};
% f     : a 2X1 cell array with body forcing functions 
%                                  f{1} = fx = @(x,y,t) ...;  for u - eqt
%                                  f{2} = fy = @(x,y,t) ...;  for v - eqt
% flags : a cell array with flags = {PLOT, GRID, incr};
% Soln  : a cell array with the manufactured solution Soln = {u, v, p};

% unpack flags
PLOT = flags{1};
GRID = flags{2};
incr = flags{3};
SAVE = flags{4};
tol = flags{5};

% unpack physical parameters and characteristic scales
%prms = [mu, rho, Re, Uchar, Xchar, Tchar, Prchar, Bchar, ...
%        Plchar, Dp, Pe, Rchar,Pmax];
mu = prms(1);
rho = prms(2);
Re = prms(3);
Uchar = prms(4);
Xchar = prms(5);
Tchar = prms(6);
Prchar = prms(7);
Bchar = prms(8);
Plchar = prms(9);
Dp = prms(10);
Pe = prms(11);
Rchar = prms(12);
Pb_max = prms(13);

% Unpack spatial parameters 
Nx = D(1);     
Ny = D(2);
D(3:8) = D(3:8)/Xchar; % Nondimensionalize lengths of H
WL = D(3); Wmid = D(4); WR = D(5); 
Hlow  = D(6); Hmid = D(7); Hup = D(8);

% Get grid for FVM 
[XU,YU,XV,YV,xp,yp,N,mindx,mindy] = ...
                  NSGridH(Nx,Ny,WL,Wmid,WR,Hlow,Hmid,Hup,GRID,PLOT);

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Create Mesh Grid for Plotting
[XP, YP] = meshgrid(xp(2:Nx+1),yp(2:Ny+1)); %Ignore ghost cells
[XC, YC] = meshgrid(xp,yp); % Grid for concentrations with ghost cells

% Unpack Flow Rates for IC and Caclulate Inlet Velocity Profile (BC for V)
Q1 = Q(1,1);  Q3 = Q(1,2);
BC_in = inletBC_V(Q1,Q3,WL,Wmid,WR,Xchar,Uchar,XV,N);
maxUV = max(abs(BC_in));

% Nondimensionalize inlet/outlet pressures
P_io(1:4) = P_io(1:4)/Prchar;

% Enforce CFL cfl = Uchar Tchar * dt / (Xchar*dx) 
% Nondimensionalize time and gather time params Nt, dt, Tf
Tf = Tfdim/Tchar;
dt = 0.75/max(maxUV/mindx,maxUV/mindy); % Less restrictive
t = 0;

% Get Matrices Au, Av, Ap
Au = BuildAuH(N,XU,YU,XV,YV,Re,dt);
Av = BuildAvH(N,XU,YU,XV,YV,Re,dt);
Apr = BuildApH(N,xp,yp,XU,YV);
Apr = decomposition(Apr); 

% Initial Conditions Velocity and Pressure 
U0 = 0*XU;
V0 = 0*XV;
P0 = initial_Pr(N,D,P_io,Xchar,Prchar,XP,YP);
% load('data/UVP_PbLow_N200_Grid1.mat','U','V','P'); 
% load('data/UVP_PETSc_N257_Grid1.mat','U','V','P')
% U = U/Uchar; V = V/Uchar; P = P/Prchar;
% U0 = U; V0 = V; P0 = P;


% Initial Conditions Mobile Platelets
[Pm0, Pm0Ghst] = initial_Pm(N,XC,YC);
Pm = Pm0;
PmGhst = Pm0Ghst;

% Initial Conditions Bound Platelets
[Pb0, Pb0Ghst] = initial_Pm(N,XC,YC);
Pb = Pb0;
PbGhst = Pb0Ghst;

% Initial Brinkman Term
B0 = alpha_fxn(theta_B(Plchar*Pb0,Pb_max),Bchar);
B = alpha_fxn(theta_B(Plchar*Pb,Pb_max),Bchar);

% Initial Hindered Transport Function
W = Wfxn(theta_T(Plchar*Pm, Plchar*Pb, Pb_max));
W0 = W;
Wghst = Wfxn(theta_T(Plchar*Pm0Ghst, Plchar*Pb0Ghst, Pb_max));
W0ghst = Wghst;

% Set Body Force Terms fx and fy 
fx = @(x,y,t)  0*x.*y.*t;
fy = @(x,y,t)  0*x.*y.*t;

% k_adhesion for reactions 
% to non-dim divide by Rchar
% k_adh = k_adh_fxn(XP,YP,Xchar) * Plchar / Rchar;
k_adh = k_adh_fxn(XP,YP,Xchar) *Uchar./Plchar./Xchar;
% k_adh = zeros(Ny,Nx);
% act_rate = (2e-10/6.022); % standard slow rate
% act_rate = (2^2/6.022); % faster rate for test
% k_adh(NyL+2:NyM,NxL+2:NxM-1) = act_rate*Xchar*Plchar./Uchar; % cm^3 / s
% k_adh(NyM,NxL+2:NxM-1) = act_rate*Xchar*Plchar./Uchar; % cm^3 / s
% k_adh(NyL+2:NyM,NxL+2:NxM-1) = act_rate*Uchar/Xchar/Plchar; % 
% k_adh(NyM,NxL+2:NxM-1) = act_rate*Uchar/Xchar/Plchar; % 


% --- One step of Forward Euler ---
[NLU,NLV] = NonLin(U0,V0,N,XU,YU,XV,YV);
NLU = 2/3*NLU;
NLV = 2/3*NLV;
NLU0 = 0*NLU;
NLV0 = 0*NLV;
[WU0,WV0] = hinderedUV_UpW(N,U0,V0,W0);

[advPm] = advectC_HiRes(Pm0,WU0,WV0,N,XU,YV,Pm0Ghst,dt);

% Reactions for ADR
RR_m = 2/3 * React_Pm(Pm(2:Ny+1,2:Nx+1),Pb(2:Ny+1,2:Nx+1),k_adh,Pb_max,Rchar);
RR0_m = 0*RR_m;
RR_b = 2/3 * React_Pb(Pm(2:Ny+1,2:Nx+1),Pb(2:Ny+1,2:Nx+1),k_adh,Pb_max,Rchar);
RR0_b = 0*RR_b;
    

% Forcing Terms at t = t0
Quo = fx(XU,YU,t); 
Qvo = fy(XV,YV,t); 

% load('data/longRun_normalK_adh_t8904');
% Apr = BuildApH(N,xp,yp,XU,YV);
% Apr = decomposition(Apr); 

% Initiate Plots (during time stepping)
if PLOT == 2
    figure(1)
% Initiate Stream Line plot
elseif PLOT >2 
    figure(1)
    plt_save = [0.01027,0.0134,0.25,0.5,1,2,4,6,8,9];
    if (SAVE==1)
        plt_ctr = 1;
        vid = VideoWriter(['figs/streamLinesTf',num2str(Tchar*Tf),'.mp4'],'MPEG-4');
        open(vid);
    end
end

steady = 1;
tic
k = 1;
MatBuilds = 0;



% Step in time
while (t < Tf) %|| (steady > tol)
    k = k+1;
    t = t + dt;
    
    % Forcing Terms
    Qu = fx(XU,YU,t); 
    Qv = fy(XV,YV,t); 
    
    % Get Diagonal Matrices for Brinkman Terms Bu and Bv
    B = alpha_fxn(theta_B(Plchar*Pb,Pb_max),Bchar);
    Bu = BuildBuH(N,XU,XV,YV,B);
    Bv = BuildBvH(N,XU,YU,YV,B);  
    
    % Prediction - Solve u-momentum equation
    Fu = RHS_u_H(U0,NLU,NLU0,P0,Qu,Quo,N,dt,Re,XU,YU,XV,YV,B0);
    TMPU = (Au+Bu) \ Fu; 

    % Prediction - Solve v- momentum equation
    Fv = RHS_v_H(V0,NLV,NLV0,P0,Qv,Qvo,N,dt,Re,XU,YU,XV,YV,BC_in,B0);
    TMPV = (Av+Bv) \ Fv;

    % Reshape TMPU and TMPV into Matrices
    [U,V] = reshapeUV_H(NxL,NxM,Nx,NyL,NyM,Ny,TMPU,TMPV);
    
    % Correction - Solve PHI Equation
    Fp = RHS_p_H(N,U,V,dt,XU,YV,yp,P0,P_io);
    TMPP = Apr\Fp;

    % Reshape TMPP into Matrix
    PHI = reshapePHI_H(N,TMPP);
    
    % Correction - Get Current U, V and P
    [U,V,P] = correction(U,V,P0,PHI,N,dt,xp,yp,BC_in);
    
    % Solve for Mobile Platelets
    BC_PmIn = inletBC_Pm(Plchar,WL,Wmid,Xchar,XC,N,t,Tchar);
    Apm = BuildAcH_Hindered(N,XC,YC,XU,YV,Pe,dt,U,V,W,Wghst);
    Fpm = RHS_c_H_ADRW(N,Pm0,advPm,dt,XU,YV,XC,YC,Pe,BC_PmIn,...
                       RR_m,RR0_m,Pm0Ghst,U,V,W0,W0ghst);
    TMPC = Apm \ Fpm; 
    
    % Reshape C and get ghost nodes for next NL term and RHS vector
    Pm = reshapeC(NxL,NxM,Nx,NyL,NyM,Ny,TMPC);
    [Pm,PmGhst] = ghostNodesC(N,Pm,BC_PmIn);

    % Solve for Bound Platelets
    [Pb,PbGhst] = react_ode_solve(N,RR_b,RR0_b,Pb,dt);
%     Pb(abs(Pb)<1e-99)=0;
    
    % Check if steady-state has been reached
    if k > 50  && ~mod(k,incr)
        Du = max(max(abs(Uchar*U-Uchar*U0)));
        Dv = max(max(abs(Uchar*V-Uchar*V0)));
        steady = max([Du,Dv]);
        DIV = median( median( abs(div(U,V,XU,YV, Nx, Ny) ))); % check div(u) = 0
        steadyPm = max(max(abs(Pm - Pm0)));
        disp(['SteadyUV = ',num2str(steady),'  Div(U) = ',num2str(DIV),...
              '  SteadyPm = ',num2str(steadyPm),...
              '   t = ',num2str(Tchar*t)])
    end
    
    % Updates (fluid)
    U0 = U;
    V0 = V;
    P0 = P;
    B0 = B;
    Quo = Qu;
    Qvo = Qv;
    NLU0 = NLU;
    NLV0 = NLV;
    [NLU,NLV] = NonLin(U,V,N,XU,YU,XV,YV);
    
    % Updates Hindered Transport Function and Hindered Velocity
    W0 = Wfxn(theta_T(Plchar*Pm, Plchar*Pb, Pb_max));
    W0ghst = Wfxn(theta_T(Plchar*PmGhst, Plchar*PbGhst, Pb_max));
    W = Wfxn(theta_T(Plchar*(2*Pm-Pm0), Plchar*(2*Pb-Pb0), Pb_max));
    Wghst = Wfxn(theta_T(Plchar*(2*PmGhst-Pm0Ghst),...
                         Plchar*(2*PbGhst-Pb0Ghst), Pb_max));
%     if k == 130
%         disp(Tchar*t)
%     end
    % Update Platelets
    Pm0 = Pm;
    Pm0Ghst = PmGhst;
    Pb0 = Pb;
    Pb0Ghst = PbGhst;
    
    [WU0,WV0] = hinderedUV_UpW(N,U0,V0,W0);
    [advPm] = advectC_HiRes(Pm0,WU0,WV0,N,XU,YV,Pm0Ghst,dt);
    RR0_m = RR_m;
    RR_m = React_Pm(Pm(2:Ny+1,2:Nx+1),Pb(2:Ny+1,2:Nx+1),k_adh,Pb_max/Plchar,1);
    RR0_b = RR_b;
    RR_b = React_Pb(Pm(2:Ny+1,2:Nx+1),Pb(2:Ny+1,2:Nx+1),k_adh,Pb_max/Plchar,1);

    % Check CFL
    maxU = max(max(abs(U))); maxV = max(max(abs(V)));
    cfl = 1 / max(maxU/mindx, maxV/mindy); % Less Restrictive
    if dt > 0.8 * cfl
        MatBuilds = MatBuilds + 1;
        disp(['Rebuilding Matrices, old dt = ',num2str(dt),...
         ', new dt = ',num2str(0.75*cfl)])
        disp(['Number of rebuilds = ',num2str(MatBuilds)])
        dt = 0.75 * cfl;
        Au = BuildAuH(N,XU,YU,XV,YV,Re,dt);
        Av = BuildAvH(N,XU,YU,XV,YV,Re,dt); 
    end

     % Save figures at specified times
    if (Tchar*t<=(plt_save(plt_ctr)+Tchar*0.6*dt)) ...
            && (Tchar*t>=(plt_save(plt_ctr)-Tchar*0.6*dt))
        clf
    makeVectorfield_PbPm(Xchar*D,N,Uchar*U,Uchar*V,Xchar*XU,...
                     Xchar*YU,Xchar*XV,Xchar*YV,Xchar*XP,Xchar*YP,...
                     Plchar*Pb(2:Ny+1,2:Nx+1),Plchar*Pm(2:Ny+1,2:Nx+1),Tchar*t,max(max(Pb*Plchar)))
        saveas(gcf, ['figs/PmPmPlot',num2str(plt_ctr),'.png'])
        saveas(gcf, ['figs/PmPmPlot',num2str(plt_ctr),'.fig'])
        plt_ctr = plt_ctr + 1;
        pause(0.1)
    end

    % plot solution at various values of time
    if PLOT == 3 && ~mod(k,incr) % Plot Stream Function
        clf
        makeVectorfield_PbPm(Xchar*D,N,Uchar*U,Uchar*V,Xchar*XU,...
                         Xchar*YU,Xchar*XV,Xchar*YV,Xchar*XP,Xchar*YP,...
                         Plchar*Pb(2:Ny+1,2:Nx+1),Plchar*Pm(2:Ny+1,2:Nx+1),Tchar*t,max(max(Pb*Plchar)))
%         makeVectorfieldBrinkman(Xchar*D,N,Uchar*U,Uchar*V,Xchar*XU,Xchar*YU,...
%                                 Xchar*XV,Xchar*YV,Xchar*XP,Xchar*YP,...
%                                 Plchar*Pm(2:Ny+1,2:Nx+1),Plchar,t*Tchar)
%         makePcolorsPlatelets(Xchar*D,Xchar*XP,Xchar*YP,Tchar*t,...
%                             Plchar*Pm(2:Ny+1,2:Nx+1),Plchar,...
%                             Plchar*Pb(2:Ny+1,2:Nx+1),Pb_max)
       
        pause(0.1)
        if (SAVE==1)
            figure(1)
            frame = getframe(gcf);
            writeVideo(vid,frame);
        end
    end
    if steady < tol
        disp(['steady state reached at t = ', num2str(Tchar*t)])
        break
    end
end
CPU_time = toc

% Save a movie of Stream Lines
if (PLOT == 3) && (SAVE == 1)
    % Output the movie as an avi file
    close(vid);
end


% Return dimensional U,V,P and associated grids
S = {Uchar*U, Uchar*V, Prchar*P, Plchar*Pm, Plchar*Pb};
B = B * Bchar;
grd = {Xchar*XU,Xchar*YU,Xchar*XV,Xchar*YV,...
      Xchar*XP,Xchar*YP, Xchar*XC, Xchar*YC};
Tf = Tchar*t;
