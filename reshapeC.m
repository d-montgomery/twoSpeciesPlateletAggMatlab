function [C] = reshapeC(NxL,NxM,Nx,NyL,NyM,Ny,TMPC)
% Inputs:
% NxL,NxM,Nx,NyL,NyM,Ny are various numbers of cells in regions of the H
% TMPC is a long vector from solving A TMPC = F 

% --- Reshape C ---
C = zeros(Ny,Nx);

% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL) + (Nx - NxM); 

% Lower Legs
j = [1:NxL, NxM+1:Nx];
for i = 1:NyL
    j_ctr = (1 : Nsml) + (i-1)*Nsml;
    C(i,j) = TMPC( j_ctr ); 
end

% Middle of H
j = 1:Nx;
for i = NyL+1 : NyM
    j_ctr = NyL*Nsml + j + (i-NyL-1)*(Nx);
    C(i,j) = TMPC( j_ctr ); 
end

% Upper Legs
j = [1:NxL, NxM+1:Nx];
for i = NyM+1 : Ny
    j_ctr = NyL*Nsml + (NyM - NyL)*(Nx) + (1 : Nsml) + (i-NyM-1)*Nsml;
    C(i,j) = TMPC( j_ctr ); 
end
