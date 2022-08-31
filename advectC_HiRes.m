function[advC] = advectC_HiRes(C0,U,V,N,XU,YV,C0Ghst,dt)
% Calculates the advection terms in advection-diffusion equation
% Inputs:
% C0 = concentration from prev. time step
% U = x comp. of velocity field
% V = y comp. of velocity field
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% XU and YU is where u(x,y,t) is stored on grid
% xc and yc is where c(x,y,t) is stored on grid
% t = t^{n+1}
% BC = cell array with vectors containing boundary conditions

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Compute advection terms for c equation -----------------------------
advC = zeros(Ny,Nx);
F = zeros(Ny+2,Nx+1);
G = zeros(Ny+1,Nx+2); 

% Wash Channel
for i = 1:Ny
for j = 1:NxL
    
    % Get all C values (takes into acount ghost cells)
    CE = C0(i+1,j+2); CN = C0(i+2,j+1); CW = C0(i+1,j); CS = C0(i,j+1);
    CP = C0(i+1,j+1);
    
    if j == 1
        CWW = CW;
    else
        CWW = C0(i+1,j-1);
    end
    if i == 1
        CSS = CS; 
    else 
        CSS = C0(i-1,j+1);
    end
    
    % Calculate DCU terms
    advC(i,j) = calc_DCU(i,j,XU,YV,CE,CN,CW,CS,CP,U,V);
    
    % Calculate CTU Updates
    [F,G] = calc_CTU(i,j,XU,YV,CP,CW,CS,U,V,F,G,dt);
    
    % Calculate Second Order Limited Corrections
    [F,G] = calc_SOC(i,j,XU,YV,CE,CN,CW,CS,CP,CWW,CSS,U,V,F,G,dt);
    
end
end

% Injury Channel
for i = NyL+1:NyM
for j = NxL+1:NxM
    
    % Get all C values (takes into acount ghost cells)
    CE = C0(i+1,j+2); CN = C0(i+2,j+1); CW = C0(i+1,j); CS = C0(i,j+1);
    CP = C0(i+1,j+1);
    
    %\Gamma 6
    if (i == NyL+1) && (j == NxL+1) % Corner 2
        CS = C0Ghst(2);
    elseif (i == NyL+1) && (j == NxM) % Corner 4
        CS = C0Ghst(4);
    end
    
    %\Gamma 5
    if (i == NyM) && (j == NxL+1) % Corner 1
        CN = C0Ghst(1);
    elseif (i == NyM) && (j == NxM) % Corner 3
        CN = C0Ghst(3);
    end
    
    % Extra Nodes for calculating r in flux limiter
    CWW = C0(i+1,j-1);
    %\Gamma 6
    if (i == NyL+2) && (j == NxL+1) % Corner 2
        CSS = C0Ghst(2);
    elseif (i == NyL+2) && (j == NxM) % Corner 4
        CSS = C0Ghst(4);
    else 
        CSS = C0(i-1,j+1);
    end
    
    
    % Calculate DCU terms
    advC(i,j) = calc_DCU(i,j,XU,YV,CE,CN,CW,CS,CP,U,V);
    
    % Calculate CTU Updates
    [F,G] = calc_CTU(i,j,XU,YV,CP,CW,CS,U,V,F,G,dt);
    
    % Calculate Second Order Limited Corrections
    [F,G] = calc_SOC(i,j,XU,YV,CE,CN,CW,CS,CP,CWW,CSS,U,V,F,G,dt);
end
end

% Blood Channel
for i = 1:Ny
for j = NxM+1:Nx
    % Get all C values (takes into acount ghost cells)
    CE = C0(i+1,j+2); CN = C0(i+2,j+1); CW = C0(i+1,j); CS = C0(i,j+1);
    CP = C0(i+1,j+1);
    
    CWW = C0(i+1,j-1);
    if i == 1
        CSS = CS; 
    else 
        CSS = C0(i-1,j+1);
    end
    
    % Calculate DCU terms
    advC(i,j) = calc_DCU(i,j,XU,YV,CE,CN,CW,CS,CP,U,V);
    
    % Calculate CTU Updates
    [F,G] = calc_CTU(i,j,XU,YV,CP,CW,CS,U,V,F,G,dt);
    
    % Calculate Second Order Limited Corrections
    [F,G] = calc_SOC(i,j,XU,YV,CE,CN,CW,CS,CP,CWW,CSS,U,V,F,G,dt);
end
end



% Update advC with CTU and SOC corrections
advC(1:Ny,1:Nx) = advC(1:Ny,1:Nx)...
             + (YV(2:Ny+1,2:Nx+1) - YV(1:Ny,2:Nx+1))...
             .*(F(2:Ny+1,2:Nx+1) - F(2:Ny+1,1:Nx)) ...
             + (XU(2:Ny+1,2:Nx+1) - XU(2:Ny+1,1:Nx))...
             .*(G(2:Ny+1,2:Nx+1) - G(1:Ny,2:Nx+1));




end

function [DCU] = calc_DCU(i,j,XU,YV,CE,CN,CW,CS,CP,U,V)

% Width and height of cell
dx = XU(i+1,j+1) - XU(i+1,j);
dy = YV(i+1,j+1) - YV(i,j+1);

% East
Ue = U(i+1,j+1);
UeCe = max(0,Ue)*CP + min(0,Ue)*CE;

% West
Uw = U(i+1,j);
UwCw = max(0,Uw)*CW + min(0,Uw)*CP;


% North
Vn = V(i+1,j+1);
VnCn = max(0,Vn)*CP + min(0,Vn)*CN;

% South
Vs = V(i,j+1);
VsCs = max(0,Vs)*CS + min(0,Vs)*CP;


DCU = (UeCe - UwCw)*dy + (VnCn - VsCs)*dx;
    
end

function [F,G] = calc_CTU(i,j,XU,YV,CP,CW,CS,U,V,F,G,dt)

% dt divided by width and/or height of cell
dtdx = dt/(XU(i+1,j+1) - XU(i+1,j));
dtdy = dt/(YV(i+1,j+1) - YV(i,j+1));

% Calculate Wave Jumps
dcw = CP - CW;
dcs = CP - CS;

% Velocities
Uw_m = min(0,U(i+1,j));
Uw_p = max(0,U(i+1,j));
Vs_m = min(0,V(i,j+1));
Vs_p = max(0,V(i,j+1));

% --- F Updates
% Fsw
F(i,j) = F(i,j) - 0.5 * dtdy * Uw_m * Vs_m * dcs;

% Fse 
F(i,j+1) = F(i,j+1) - 0.5 * dtdy * Uw_p * Vs_m * dcs;

% Fw 
F(i+1,j) = F(i+1,j) - 0.5 * dtdy * Uw_m * Vs_p * dcs;

% Fe
F(i+1,j+1) = F(i+1,j+1) - 0.5 * dtdy * Uw_p * Vs_p * dcs;

% --- G Updates
% Gsw 
G(i,j) = G(i,j) - 0.5 * dtdx * Uw_m * Vs_m * dcw;

% Gnw
G(i+1,j) = G(i+1,j) - 0.5 * dtdx * Uw_m * Vs_p * dcw;

% Gs
G(i,j+1) = G(i,j+1) - 0.5 * dtdx * Uw_p * Vs_m * dcw;

% Gn 
G(i+1,j+1) = G(i+1,j+1) - 0.5 * dtdx * Uw_p * Vs_p * dcw;
end

function [F,G] = calc_SOC(i,j,XU,YV,CE,CN,CW,CS,CP,CWW,CSS,U,V,F,G,dt)
    Uw = U(i+1,j);
    Vs = V(i,j+1);
    
    % dt divided by width and/or height of cell
    dx = XU(i+1,j+1) - XU(i+1,j);
    dy = YV(i+1,j+1) - YV(i,j+1);
    dtbydx = dt/dx;
    dtbydy = dt/dy;
    
    % Calculate Wave Jumps
    dcw = CP - CW;
    dcs = CP - CS;

    % --- F SOC update
    r = 1;
    if (Uw > 0)  
        r = (CW - CWW)/dcw;
    elseif Uw < 0 
        r = (CE - CP)/dcw;
    end
    Sw = 0.5 * abs(Uw) * (1 - dtbydx*abs(Uw))*phi(r)*dcw;
    F(i+1,j) = F(i+1,j) + Sw;
    
    % --- G SOC update
    r = 1;
    if Vs > 0 
        r = (CS - CSS)/dcs;
    elseif Vs < 0 
        r = (CN - CP)/dcs;
    end

    Ss = 0.5 * abs(Vs) * (1 - dtbydy*abs(Vs))*phi(r)*dcs;
    G(i,j+1) = G(i+1,j) + Ss;
    
    % --- Second Order Transverse Wave Corrections
    F(i+1,j+1) = F(i+1,j+1) + dtbydy * max(Uw,0) * Ss; % Fe
    F(i+1,j) = F(i+1,j) + dtbydy * min(Uw,0) * Ss; % Fw
    F(i,j+1) = F(i,j+1) - dtbydy * max(Uw,0) * Ss; % Fse
    F(i,j) = F(i,j) - dtbydy * min(Uw,0) * Ss; % Fsw
    
    G(i+1,j+1) = G(i+1,j+1) + dtbydx * max(Vs,0) * Sw; %Gn
    G(i,j+1) = G(i,j+1) + dtbydx * min(Vs,0) * Sw; % Gs
    G(i+1,j) = G(i+1,j) - dtbydx * max(Vs,0) * Sw; %Gnw
    G(i,j) = G(i,j) - dtbydx * min(Vs,0) * Sw; %Gsw
end



function lim = phi(r)

% % MC Limiter
% lim = max(0, min([(1+r)/2, 2, 2*r]));

% Minmod
lim = max(0,min(1,r));

end 