function[Fp] = RHS_p_H(N,U,V,dt,XU,YV,yp,P0,P_io)
% Inputs:
% N = [NxL NxM Nx NyL NyM Ny] for H domain
% U = predictive velocity u*
% V = predictive velocity v*
% dt = temporal step size
% XU = x-comp of mesh for U
% YV = y-comp of mesh for V
% yp = y-comp of mesh for P
% P0 = p^n pressure from prev time step
% P_io = [Pw_in Pw_out Pb_in Pw_out P_fix] inlet/outlet pressures
%        P_fix = 1 or 0:  1 -> fixes pressure at outlet

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Get pressures at the inlet/outlet
Pw_in = P_io(1);
Pw_out = P_io(2);
Pb_in = P_io(3);
Pb_out = P_io(4);
P_fix = P_io(5); % This flag determines if dphi/dy = 0 or P_fix - P^n 

% --- Matrix Sizes ---
% Ap is NNXY x NNXY matrix without augmentation
NNXY = (Nx)*(Ny) - (NxM - NxL) * ( Ny - NyM + NyL);% Size of sml sub matrices = Nsml x Nsml
Nsml = (NxL) + (Nx - NxM); 

Fp = zeros(NNXY,1);

% ---- Fill Cells in Lower Legs of H ------------------------------
for i = 1:NyL
    j_ctr = NxL;
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*(i-1)+j_ctr;
    end
    
    BC2 = 0;
    BC4 = 0;
    dx = XU(i+1,j+1) - XU(i,j);  
    dy = YV(i+1,j+1) - YV(i,j+1);
    
    if (P_fix == 1) && (i == 1) 
        as = 1/(yp(i+1)-yp(i))/dy; 
        if (j < NxL+1 )
            Pw_n = (3*P0(i,j) - P0(i+1,j))/2;
            BC2 = -2*as*( Pw_out - Pw_n);
        else
              Pb_n = (3*P0(i,j) - P0(i+1,j))/2;
              BC4 = -2*as*( Pb_out - Pb_n );
        end
    end
   
    Fp( ROW_ind ) = 1/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx ...
                          +  ( V(i+1,j+1) - V(i,j+1) )/dy ) + BC2 + BC4;
    
end  
end

% ---- Fill Cells in Injury Channel of H ------------------------------
for i = NyL+1 : NyM   
for j = 1:Nx
    ROW_ind = Nsml*(NyL) + (i - NyL -1) * Nx + j;
    
    dx = XU(i+1,j+1) - XU(i,j);  
    dy = YV(i+1,j+1) - YV(i,j+1);
   
    Fp( ROW_ind ) = 1/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx ...
                          +  ( V(i+1,j+1) - V(i,j+1) )/dy );
    
end  
end


% ---- Fill Cells in Upper Legs of H ------------------------------
for i = NyM+1:Ny
    j_ctr = NxL;
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j_ctr;
    end
    
    dx = XU(i+1,j+1) - XU(i,j);  
    dy = YV(i+1,j+1) - YV(i,j+1);
    
    BC1 = 0;
    BC3 = 0;
    
%     if (P_fix == 1) && (i == Ny) 
%         an = 1/(yp(i+2)-yp(i+1))/dy; 
%         if (j < NxL+1 )
% %             Pw_n = P0(i-1,j) + (YV(Ny+1,j) - yp(i))/(yp(i+1) - yp(i))...
% %                      *(P0(i,j) - P0(i-1,j) );
%             Pw_n = (3*P0(i-1,j) - P0(i,j))/2;
%             BC1 = -2*an*( Pw_in - Pw_n);
%         else
% %             Pb_n = P0(i-1,j) + (YV(Ny+1,j) - yp(i))/(yp(i+1) - yp(i))...
% %                      *(P0(i,j) - P0(i-1,j) );
%             Pb_n = (3*P0(i-1,j) - P0(i,j))/2;
%             BC3 = -2*an*( Pb_in - Pb_n );
%         end
%     end
   
    Fp( ROW_ind ) = 1/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx ...
                          +  ( V(i+1,j+1) - V(i,j+1) )/dy ) + BC1 + BC3;
    
end  
end

% % Augmentation
% Fp(NNXY+1) = 0;





