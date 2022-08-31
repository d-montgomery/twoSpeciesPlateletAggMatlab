function [PHI] = reshapePHI_H(N,TMPP)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Reshape PHI ---
PHI = zeros(Ny,Nx);

% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL) + (Nx - NxM); 

% Lower Legs
j = [1:NxL, NxM+1:Nx];
for i = 1:NyL
    j_ctr = (1 : Nsml) + (i-1)*Nsml;
    PHI(i,j) = TMPP( j_ctr ); 
end

% Middle of H
j = 1:Nx;
for i = NyL+1 : NyM
    j_ctr = NyL*Nsml + j + (i-NyL-1)*(Nx);
    PHI(i,j) = TMPP( j_ctr ); 
end

% Upper Legs
j = [1:NxL, NxM+1:Nx];
for i = NyM+1 : Ny
    j_ctr = NyL*Nsml + (NyM - NyL)*(Nx) + (1 : Nsml) + (i-NyM-1)*Nsml;
    PHI(i,j) = TMPP( j_ctr ); 
end
