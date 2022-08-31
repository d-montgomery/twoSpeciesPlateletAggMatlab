function[F] = RHS_v_H(V0,NLV,NLV0,P0,Qv,Qvo,N,dt,Re,XU,YU,XV,YV,BC_in,B)
% INPUTS:
% U0 = Previous velocity U
% NLU = Nonlinear terms at t = t_n
% NLU0 = Nonlinear terms at t = t_{n-1}
% P0 = Pressure at t = t_n
% Qv = Body forcing terms at t = t_{n+1}
% Qvo = Body forcing terms at t = t_n
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% dt = temporal step size
% Re = Reynolds number
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% cd are the coefficients for the inlet/outlet BC's
% BC = cell array with vectors containing boundary conditions
% t = t^{n+1}


% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Av is NNXY x NNXY matrix
NNXY = (Nx+2)*(Ny+1) - (NxM - NxL -2) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 2) + (Nx + 2 - NxM); 

F = zeros( NNXY,1 );

% --- Boundary Conditions ------------------------------------------------

% Top Left BC (\partial \Omega_1)
for j = 1:NxL+2
    F(NNXY-Nsml+j) = BC_in(j); 
end

% Top Right BC (\partial \Omega_3)
j_ctr = NxL+2; % indexing for right bc starts at NxL + 4
for j = NxM+1 : Nx+2
    j_ctr = j_ctr + 1;
    F(NNXY-Nsml+j_ctr) = BC_in(j);
end

% % Bottom Boundaries (\partial \Omega_i, i = 2, 4)
% % c_i * Vp + d_i * (dVp/dt + Vinf * dVp/dy) = h_i
% j_ctr = NxL+2; % indexing for right bc starts at NxL + 4
% for j = [1 : NxL+2, NxM+1 : Nx+2]
%     if j < NxL + 3 % bottom left bc
%         F(j) = V0(1,j) / dt;
%      
%     else % bottom right bc
%         % Principal
%         j_ctr = j_ctr +1;
%         F(j_ctr) = V0(1,j) / dt;
%     end
% end

% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if j < NxM+2
        ROW = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    end
    
    val = vRHS_val(i,j,XU,YU,XV,YV,V0,P0,B,NLV,NLV0,Qv,Qvo,Re,dt);

    F( ROW ) = val;
  
end  
end


% ---- Fill Inner cells at Lower Interface yv_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if (i == NyL) && (j > NxM+1)
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    else
        ROW = Nsml*(i-1)+j;
    end
    
    val = vRHS_val(i,j,XU,YU,XV,YV,V0,P0,B,NLV,NLV0,Qv,Qvo,Re,dt);

    F( ROW ) = val;
  
end  
end


% ---- Fill inner cells in middle of H for i = NyL +2, ..., NyM
for i = NyL+2 : NyM
for j = 2:Nx+1
    
    ROW = NyL*Nsml + (i - NyL-1)*(Nx+2) + j;
    
    val = vRHS_val(i,j,XU,YU,XV,YV,V0,P0,B,NLV,NLV0,Qv,Qvo,Re,dt);

    F( ROW ) = val;
  
end  
end


% ---- Fill Inner cells at Upper Interface yv_i, i = NyM+1, NyM+2 ----
for i = NyM+1 : NyM+2
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if (i == NyM+2) && (j > NxM+1)
        j_ctr = j_ctr + 1;
        ROW = Nsml*NyL + (NyM+1 - NyL) * (Nx+2) + j_ctr;
    else
        ROW = Nsml*NyL + (i-1 - NyL) * (Nx+2)+j;
    end
    
    val = vRHS_val(i,j,XU,YU,XV,YV,V0,P0,B,NLV,NLV0,Qv,Qvo,Re,dt);

    F( ROW ) = val;
  
end  
end


% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+3 : Ny
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if j < NxM+2
        ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j_ctr;
    end
    
    val = vRHS_val(i,j,XU,YU,XV,YV,V0,P0,B,NLV,NLV0,Qv,Qvo,Re,dt);

    F( ROW ) = val;
  
end  
end

end

function val = vRHS_val(i,j,XU,YU,XV,YV,V0,P0,B,NLV,NLV0,Qv,Qvo,Re,dt)
    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);
    dy = YU(i+1,j)-YU(i,j);
    
    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);
    
    % Velocities
    VP = V0(i,j);
    VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);
   
    % Diffusion Terms
    dw = -1/Re * ( VP - VW )*dy/(XP-XW)/2;
    de = 1/Re * ( VE - VP )*dy/(XE-XP)/2;
    ds = -1/Re * ( VP - VS )*dx/(YP-YS)/2;
    dn = 1/Re * ( VN - VP )*dx/(YN-YP)/2;
    D = dw + de + ds + dn;

    % Pressure
    Pn = P0(i,j-1); Ps = P0(i-1,j-1);
    PR = (Pn-Ps) * dx;

    % Brinkman Term 
    [Bp,~,~] = Bp_Vcoeff(i,j,XU,YU,YV,B);

    val = (1/dt - Bp/2)*V0(i,j)*dx*dy - 1.5*NLV(i,j) + ...
                              + 0.5*NLV0(i,j) + D - PR + ...
                              + 0.5*(Qv(i,j) + Qvo(i,j))*dy*dx;
end