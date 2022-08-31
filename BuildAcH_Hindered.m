function[A] = BuildAcH_Hindered(N,XC,YC,XU,YV,Pe,dt,U,V,W,Wghst)
% Builds the pentadiagonal matrix for advection-diffusion equation
% Inputs:
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% xc and yc is where c(x,y,t) is stored on grid
% XU and YV are the locations of the cell walls
% Pe is the Peclet number
% dt is the temporal time step size
% W is the hindered transport function evaluated at 2*C^{n} - C0^{n-1}
% Wghst contains the values of W for the Gamma 5 and 6 boundaries 

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Ac is NNXY x NNXY matrix 
NNXY = (Nx)*(Ny) - (NxM - NxL) * ( Ny - NyM + NyL);

% Size of sml sub matrices = Nsml x Nsml
Nsml = (NxL) + (Nx - NxM); 

% Initialize the counter used in building the sparse matrices
ctr = 0;

% ---- Fill Cells in Lower Legs of H ------------------------------
for i = 1:NyL
    
    j_ctr = NxL; % Starting value for index when j > NxL
    
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*(i-1)+j_ctr;
    end

    % Get Coefficients AE,AN,AW,AS,AP
    [AE,AN,AW,AS,AP] = coeffs_Ac(i,j,XU,YV,XC,YC,dt,Pe,U,V,W);
    
    % Extra terms for Ghost Nodes
    Rp = 0; Lp = 0; Bp = 0; 
    if (i == 1) % On \partial \Omega_2 or 4 boundary
        Bp = -AS;
    end
    if (j == NxL) && (i < NyL+1)  % On \Gamma_2 boundary
        Rp =  -AE;
    end
    if (j == NxM+1) && (i < NyL+1) % On \Gamma_4 boundary
        Lp = - AW;
    end 
    if (j == 1) % On \Gamma_7 boundary
        Lp = - AW;
    end 
    if (j == Nx) % On \Gamma_8 boundary
        Rp = -AE;
    end
    
    AP = AP + Lp + Rp + Bp;
    
    % Principal
    ctr = ctr+1;
    row(ctr) = ROW_ind;
    col(ctr) = row(ctr);
    A(ctr)   = AP;
    
    % East - Wash Channel
    if j < NxL
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = -AE;
    end
    
    % East - Blood Channel
    if (j > NxM) && (j < Nx)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = -AE;
    end

    % West - Wash Channel
    if (j > 1) && (j < NxL+1)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = -AW;
    end
    
    % West - Blood Channel
    if (j > NxM+1) 
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = -AW;
    end
    
    % North
    if (i == NyL) && (j > NxL) % At interface of lower legs & injury chan.
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nx;
        A(ctr)   = -AN;
    else
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nsml;
        A(ctr)   = -AN;
    end
    
    % South
    if i > 1
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nsml;
        A(ctr)   = -AS;
    end
   
end
end

% ---- Fill Cells in Injury Channel ---
for i = NyL+1 : NyM   
for j = 1:Nx
    ROW_ind = Nsml*(NyL) + (i - NyL -1) * Nx + j;

    % Width and height of cell
    dx = XU(i+1,j+1) - XU(i+1,j);
    dy = YV(i+1,j+1) - YV(i,j+1);
    
    % Location of Nodes (note xc and yc have ghost nodes)
    XE = XC(i+1,j+2);  XP = XC(i+1,j+1);  XW = XC(i+1,j);
    YN = YC(i+2,j+1);  YP = YC(i+1,j+1);  YS = YC(i,j+1);
    
    % Interpolate Hindered Transport to cell walls
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
    
    % Coefficients for Diffusion
    AE = dy*We / (2*Pe * (XE - XP));
    AN = dx*Wn / (2*Pe * (YN - YP));
    AW = dy*Ww / (2*Pe * (XP - XW));
    AS = dx*Ws / (2*Pe * (YP - YS));
    AP = dx*dy/dt + AE + AN + AW + AS;
    
    % Extra weights for ghost nodes
    Tp = 0; Bp = 0; Lp = 0; Rp = 0;
    if (i == NyM) && ((j > NxL) && (j < NxM+1)) % On \Gamma_5 boundary
        Tp = -AN;
    end 
    if (i == NyL+1) && ((j > NxL) && (j < NxM+1)) % On \Gamma_6 boundary
        Bp = - AS;
    end 
    if (j == 1) % On \Gamma_7 boundary
        Lp = - AW;
    end 
    if (j == Nx) % On \Gamma_8 boundary
        Rp = -AE;
    end
    
    AP = AP + Tp + Bp + Lp + Rp;
    
    % Principal
    ctr = ctr+1;
    row(ctr) = ROW_ind;
    col(ctr) = row(ctr);
    A(ctr)   = AP;
    
    % East 
    if j < Nx
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = -AE;
    end
    
    % West 
    if (j > 1) 
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = -AW;
    end
    
    % North
    if (i < NyM) % in injury channel
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nx;
        A(ctr)   = -AN;
    elseif (i == NyM) && (j < NxL+1) % at interface left
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nx;
        A(ctr)   = -AN;
    elseif (i == NyM) && (j > NxM) % at interface right
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nsml;
        A(ctr)   = -AN;
    end
    
    % South
    if (i > NyL+1) % in injury channel
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nx;
        A(ctr)   = -AS;
    elseif (i == NyL+1) && (j < NxL+1)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nsml;
        A(ctr)   = -AS;
    elseif (i == NyL+1) && (j > NxM)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nx;
        A(ctr)   = -AS;
    end
 
end
end

% ---- Fill Cells in Upper Legs of H ------------------------------
for i = NyM+1:Ny
    
    j_ctr = NxL; % Starting value for index when j > NxL
    
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j_ctr;
    end

    % Get Coefficients AE,AN,AW,AS,AP
    [AE,AN,AW,AS,AP] = coeffs_Ac(i,j,XU,YV,XC,YC,dt,Pe,U,V,W);
    
    % Extra weights from Ghost Nodes
    Lp = 0; Tp = 0; Rp = 0; Bp = 0;

    if (i == Ny) % On \partial \Omega_1 or 3
        Tp = AN;
    end 
    if (j == NxL) && (i > NyM)  % On \Gamma_1 boundary
        Rp =  -AE;
    end
    if (j == NxM+1) && (i > NyM) % On \Gamma_3 boundary
        Lp = - AW;
    end 
    if (j == 1) % On \Gamma_7 boundary
        Lp = - AW;
    end 
    if (j == Nx) % On \Gamma_8 boundary
        Rp = -AE;
    end

    % Update AP with Ghost Node Weights
    AP = AP + Rp + Tp + Lp + Bp;

    % Principal
    ctr = ctr+1;
    row(ctr) = ROW_ind;
    col(ctr) = row(ctr);
    A(ctr)   = AP;
    
    % East - Wash Channel
    if j < NxL
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = -AE;
    end
    
    % East - Blood Channel
    if (j > NxM) && (j < Nx)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = -AE;
    end

    % West - Wash Channel
    if (j > 1) && (j < NxL+1)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = -AW;
    end
    
    % West - Blood Channel
    if (j > NxM+1) 
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = -AW;
    end
    
    % North
    if (i < Ny)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nsml;
        A(ctr)   = -AN;
    end
    
    % South
    if (i == NyM+1) && (j < NxL+1) % at Interface of Inj and Upper Legs
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nx;
        A(ctr)   = -AS;
    else
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nsml;
        A(ctr)   = -AS;
    end
   
end
end


% Create Sparse Matrix of A
A = sparse(row,col,A,NNXY,NNXY);

end 

function [AE,AN,AW,AS,AP] = coeffs_Ac(i,j,XU,YV,XC,YC,dt,Pe,U,V,W)

% Width and height of cell
dx = XU(i+1,j+1) - XU(i+1,j);
dy = YV(i+1,j+1) - YV(i,j+1);

% Location of Nodes (note xc and yc have ghost nodes)
XE = XC(i+1,j+2);  XP = XC(i+1,j+1);  XW = XC(i+1,j);
YN = YC(i+2,j+1);  YP = YC(i+1,j+1);  YS = YC(i,j+1);

% Interpolate W to cell edges
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

% Coefficients
AE = dy*We / (2*Pe * (XE - XP));
AN = dx*Wn / (2*Pe * (YN - YP));
AW = dy*Ww / (2*Pe * (XP - XW));
AS = dx*Ws / (2*Pe * (YP - YS));
AP = dx*dy/dt + AE + AN + AW + AS;

end

function s = sgn(x)
s = -1 + heaviside(x);
end