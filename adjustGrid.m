function [U,YU,V,XV] = adjustGrid(N,U,YU,V,XV)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Fix grid to be mesh grid for U and V
U(NyL+1,NxL+2:NxM) = - U(NyL+2,NxL+2:NxM);
U(NyM+2,NxL+2:NxM) = U(NyM+1,NxL+2:NxM);
YU(NyL+1,NxL+2:NxM) = YU(NyL+1,NxL+1);
YU(NyM+2,NxL+2:NxM) = YU(NyM+2,NxL+1);

i = [1:NyL, NyM+2:Ny+1];
V(i,NxL+2) = - V(i,NxL+1);
V(i,NxM+1) = - V(i,NxM+2);
XV(i,NxL+2) = XV(NyL+1,NxL+2);
XV(i,NxM+1) = XV(NyL+1,NxM+1); 
