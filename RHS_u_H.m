function[F] = RHS_u_H(U0,NLU,NLU0,P0,Qu,Quo,N,dt,Re,XU,YU,XV,YV,B)
% INPUTS:
% U0 = Previous velocity U
% NLU = Nonlinear terms at t = t_n
% NLU0 = Nonlinear terms at t = t_{n-1}
% P0 = Pressure at t = t_n
% Qu = Body forcing terms at t = t_{n+1}
% Quo = Body forcing terms at t = t_n
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% dt = temporal step size
% Re = Reynolds number
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% ab are the coefficients for the inlet/outlet BC's
% BC = cell array with vectors containing boundary conditions
% t = t^{n+1}

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Au is NNXY x NNXY matrix, Fu is an NNXY x 1 vector
NNXY = (Nx+1)*(Ny+2) - (NxM - NxL -1) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 1) + (Nx + 1 - NxM); 

F = zeros( NNXY,1 );


% % Bottom Left Boundary (\partial \Omega_2)
% for j = 1:NxL+1
%     F(j) = U0(1,j) / dt;
% end
% 
% % Bottom Right Boundary (\partial \Omega_4)
% j_ctr = NxL+1;
% for j = NxM+1:Nx+1
%     j_ctr = j_ctr+1;
%     F(j_ctr) = U0(1,j) / dt;
% end

% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    j_ctr = NxL+2;
    
for j = [2:NxL, NxM+2:Nx]
    
    if j < NxM+2
        ROW = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    end
    
    val = uRHS_val(i,j,XU,YU,XV,YV,U0,P0,B,NLU,NLU0,Qu,Quo,Re,dt);

        F( ROW ) = val;
  
end  
end

% ---- Fill Inner cells at Lower Interface yu_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
    j_ctr = NxL+2; % Starting matrix value right H when i == NyL only 
    
    for j = [2:NxL, NxM+2:Nx]
        if (i == NyL) && j > NxM+1
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        else
            ROW = Nsml*(i-1)+j;
        end
        
        val = uRHS_val(i,j,XU,YU,XV,YV,U0,P0,B,NLU,NLU0,Qu,Quo,Re,dt);

        F( ROW ) = val;

    end 
end


% ---- Fill inner cells in middle of H for i = NyL+2, ..., NyM+1 ---
for i = NyL+2 : NyM+1
    for j = 2:Nx
        
        ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;
        
        val = uRHS_val(i,j,XU,YU,XV,YV,U0,P0,B,NLU,NLU0,Qu,Quo,Re,dt);

        F( ROW ) = val;

    end 
end

% ---- Fill Inner cells at Upper Interface yu_i, i = NyM+2, NyM+3 ----
for i =  NyM+2 : NyM+3
    j_ctr = NxL+2; % Starting matrix value right H when i == NyM+3 only 

    for j = [2:NxL, NxM+2:Nx]
        if (i == NyM+3) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + j_ctr;
        else
            ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;
        end
        
        val = uRHS_val(i,j,XU,YU,XV,YV,U0,P0,B,NLU,NLU0,Qu,Quo,Re,dt);

        F( ROW ) = val;
    end 
end

% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+4 : Ny+1
    j_ctr = NxL+2; % Starting matrix value right H when i == NyM+3 only 

    for j = [2:NxL, NxM+2:Nx]
        if j < NxM+2
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j_ctr;
        end
        
        val = uRHS_val(i,j,XU,YU,XV,YV,U0,P0,B,NLU,NLU0,Qu,Quo,Re,dt);

        F( ROW ) = val;
    end 
end

end

function val = uRHS_val(i,j,XU,YU,XV,YV,U0,P0,B,NLU,NLU0,Qu,Quo,Re,dt)
    % Cell height/width
    dx = XV(i,j+1) - XV(i,j);
    dy = YV(i,j) - YV(i-1,j);

    % Location of Nodes
    XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
    YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);

    % Velocities           
    UP = U0(i,j);
    UE = U0(i,j+1); UN = U0(i+1,j); UW = U0(i,j-1); US = U0(i-1,j);


    % Diffusion Terms
    dw = -1/Re * (UP - UW)*dy/(XP-XW)/2;
    de = 1/Re * (UE - UP)*dy/(XE-XP)/2;
    ds = -1/Re * (UP - US)*dx/(YP-YS)/2;
    dn = 1/Re * (UN - UP)*dx/(YN-YP)/2;
    D = de + dn + dw + ds;

    % Pressure Terms
    Pe = P0(i-1,j); Pw = P0(i-1,j-1);
    PR =  (Pe-Pw)* dy;

    % Brinkman Term
    [Bp,~,~] = Bp_Ucoeff(i,j,XU,XV,YV,B);

    val = (1/dt - Bp/2)*U0(i,j)*dx*dy - 1.5*NLU(i,j) +...
                                  0.5*NLU0(i,j) + D - PR + ...
                                 + (Qu(i,j) + Quo(i,j))*dx*dy/2;
end