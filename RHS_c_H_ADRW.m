function [F] = RHS_c_H_ADRW(N,C0,advC,dt,XU,YV,XC,YC,...
                                 Pe,BC,R,R0,C0Ghst,U,V,W,Wghst)
% Calculates the RHS vector in advection-diffusion-reaction equation
% Inputs:
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% C0 = concentration from prev. time step
% advC & advC0 are the advection terms from the nth and (n-1) time steps
% Q and Qo are the forcing terms from the current n+1 and nth time steps
% dt is the temporal time step size
% XU and YV are the locations of the cell walls
% xc and yc is where c(x,y,t) is stored on grid
% Pe is the Peclet number
% BC = cell array with vectors containing boundary conditions
% t = t^{n+1}
% injGhst = 4X1 vector with vals for ghost nodes at corners of \Gamma 5,6

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Ac is NNXY x NNXY matrix 
NNXY = (Nx)*(Ny) - (NxM - NxL) * ( Ny - NyM + NyL);

% Size of sml sub matrices = Nsml x Nsml
Nsml = (NxL) + (Nx - NxM); 

% Initialize RHS Vector
F = zeros(NNXY,1);

% ---- Fill Cells in Lower Legs of H ------------------------------
for i = 1:NyL
    
    j_ctr = NxL; % Starting value for index when j > NxL
    
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    end
    
    % Get all C values (takes into acount ghost cells)
    CE = C0(i+1,j+2); CN = C0(i+2,j+1); CW = C0(i+1,j); CS = C0(i,j+1);
    CP = C0(i+1,j+1);
    
    % Get Reaction values
    RP = R(i,j); RP0 = R0(i,j);
    
    % Get values values for RHS on interior cells and AE,AN,AS,AW
    [RHS] = cRHS_constD_val(i,j,CE,CN,CW,CS,CP,XU,YV,XC,YC,...
                                        dt,advC,Pe,RP,RP0,U,V,W);
    
    F( ROW ) = RHS;
    
end
end
    
% ---- Fill Cells in Injury Channel ---
for i = NyL+1 : NyM   
for j = 1:Nx
    ROW = Nsml*(NyL) + (i - NyL -1) * Nx + j;
    
     % Location of Nodes (note xc and yc have ghost nodes)
    XE = XC(i+1,j+2);  XP = XC(i+1,j+1);  XW = XC(i+1,j);
    YN = YC(i+2,j+1);  YP = YC(i+1,j+1);  YS = YC(i,j+1);
    
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

    % Interpolate Hindered Transport to cell walls
    dx = XU(i+1,j+1) - XU(i+1,j);
    dy = YV(i+1,j+1) - YV(i,j+1);
    ig = i+1; jg = j+1;
    WP = W(ig,jg); 
    We = (dx/2/(XE-XP))*W(ig,jg+1) + (XU(i+1,j+1)-XE)/(XP-XE)*WP;
    Wn = (dy/2)/(YN-YP)*W(ig+1,jg) + (YV(i+1,j+1)-YN)/(YP-YN)*WP;
    Ww = (-dx/2)/(XW-XP)*W(ig,jg-1) + (XU(i+1,j)-XW)/(XP-XW)*WP;
    Ws = (-dy/2)/(YS-YP)*W(ig-1,jg) + (YV(i,j+1)-YS)/(YP-YS)*WP;
%     We = max(0,sgn(U(i+1,j+1)))*WP + min(0,sgn(U(i+1,j+1)))*W(ig,jg+1);
%     Wn = max(0,sgn(V(i+1,j+1)))*WP + min(0,sgn(V(i+1,j+1)))*W(ig+1,jg);
%     Ww = max(0,sgn(U(i+1,j)))*W(ig,jg-1) + min(0,sgn(U(i+1,j)))*WP;
%     Ws = max(0,sgn(V(i,j+1)))*W(ig-1,jg) + min(0,sgn(V(i,j+1)))*WP;
  
    % Corners
    if (i == NyM) && (j == NxL+1) % Corner 1
        Wn = (dy/2)/(YN-YP)*Wghst(1) + (YV(i+1,j+1)-YN)/(YP-YN)*WP;
    end
    if (i == NyL+1) && (j == NxL+1) % Corner 2
        Ws = (-dy/2)/(YS-YP)*Wghst(2) + (YV(i,j+1)-YS)/(YP-YS)*WP;
    end
    if (i == NyM) && (j == NxM) % Corner 3
        Wn = (dy/2)/(YN-YP)*Wghst(3) + (YV(i+1,j+1)-YN)/(YP-YN)*WP;
    end
    if (i == NyL+1) && (j == NxM) % Corner 4
        Ws = (-dy/2)/(YS-YP)*Wghst(4) + (YV(i,j+1)-YS)/(YP-YS)*WP;
    end
    
    % Calculate Difussion Term
    D =  We*dy*(CE - CP)/(XE - XP) - Ww*dy*(CP - CW)/(XP - XW)...
        + Wn*dx*(CN - CP)/(YN - YP) - Ws*dx*(CP - CS)/(YP - YS);
    
    % Get Reaction values
    RP = R(i,j); RP0 = R0(i,j);
    
    % Reaction term
    RT = dx * dy * ( 1.5*RP - 0.5*RP0 );

    % Value for RHS vector
    RHS = CP*dx*dy/dt - advC(i,j) + 0.5*D/Pe + RT;
    
    F( ROW ) = RHS;
 
end
end


% ---- Fill Cells in Upper Legs of H ------------------------------
for i = NyM+1:Ny
    
    j_ctr = NxL; % Starting value for index when j > NxL
    
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j_ctr;
    end
    
    % Get all C values (takes into acount ghost cells)
    CE = C0(i+1,j+2); CN = C0(i+2,j+1); CW = C0(i+1,j); CS = C0(i,j+1);
    CP = C0(i+1,j+1);
    
    % Get Reaction values
    RP = R(i,j); RP0 = R0(i,j);
    
    % Get values values for RHS on interior cells and AE,AN,AS,AW
    [RHS] = cRHS_constD_val(i,j,CE,CN,CW,CS,CP,XU,YV,XC,YC,...
                                        dt,advC,Pe,RP,RP0,U,V,W);
    
    % Initialize Extra Coeffs from Ghost Nodes
    Tp = 0; 
   
    if (i == Ny) && (j < NxL + 1 ) % On \partial \Omega_1 
        dx = XU(i+1,j+1) - XU(i+1,j);
        YN = YC(i+2);  YP = YC(i+1);
        AN = dx*Wn / (2*Pe * (YN - YP));
        Tp = 2*AN*BC(j+1);
    end 
    if (i == Ny) && (j > NxM ) % On \partial \Omega_3
        dx = XU(i+1,j+1) - XU(i+1,j);
        YN = YC(i+2);  YP = YC(i+1);
        AN = dx*Wn / (2*Pe * (YN - YP));
        Tp = 2*AN*BC(j+1);
    end 
   
    
    F( ROW ) = RHS + Tp;
    
end
end

end


function [RHS]=cRHS_constD_val(i,j,CE,CN,CW,CS,CP,XU,YV,...
                                         XC,YC,dt,advC,Pe,R,R0,U,V,W)

% Width and height of cell
dx = XU(i+1,j+1) - XU(i+1,j);
dy = YV(i+1,j+1) - YV(i,j+1);

% Location of Nodes (note xc and yc have ghost nodes)
XE = XC(i+1,j+2);  XP = XC(i+1,j+1);  XW = XC(i+1,j);
YN = YC(i+2,j+1);  YP = YC(i+1,j+1);  YS = YC(i,j+1);

% Interpolate W to Cell Edges 
ig = i+1; jg = j+1;
WP = W(ig,jg); 
We = (dx/2/(XE-XP))*W(ig,jg+1) + (XU(i+1,j+1)-XE)/(XP-XE)*WP;
Wn = (dy/2)/(YN-YP)*W(ig+1,jg) + (YV(i+1,j+1)-YN)/(YP-YN)*WP;
Ww = (-dx/2)/(XW-XP)*W(ig,jg-1) + (XU(i+1,j)-XW)/(XP-XW)*WP;
Ws = (-dy/2)/(YS-YP)*W(ig-1,jg) + (YV(i,j+1)-YS)/(YP-YS)*WP;
% We = max(0,sgn(U(i+1,j+1)))*WP + min(0,sgn(U(i+1,j+1)))*W(ig,jg+1);
% Wn = max(0,sgn(V(i+1,j+1)))*WP + min(0,sgn(V(i+1,j+1)))*W(ig+1,jg);
% Ww = max(0,sgn(U(i+1,j)))*W(ig,jg-1) + min(0,sgn(U(i+1,j)))*WP;
% Ws = max(0,sgn(V(i,j+1)))*W(ig-1,jg) + min(0,sgn(V(i,j+1)))*WP;

D =  We*dy*(CE - CP)/(XE - XP) + Wn*dx*(CN - CP)/(YN - YP) ...
    -Ww*dy*(CP - CW)/(XP - XW) - Ws*dx*(CP - CS)/(YP - YS);
  

% Reaction term
RT = dx * dy * ( 1.5*R - 0.5*R0 );

% Value for RHS vector
RHS = CP*dx*dy/dt - advC(i,j) + 0.5*D/Pe + RT;

end

function s = sgn(x)
s = -1 + heaviside(x);
end