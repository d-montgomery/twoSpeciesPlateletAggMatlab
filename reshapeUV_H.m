function [U,V] = reshapeUV_H(NxL,NxM,Nx,NyL,NyM,Ny,TMPU,TMPV)



% --- Reshape U ---
U = zeros(Ny+2,Nx+1);
% Size of sml sub matrices in Au is Nsub x Nsub
Nsml = (NxL + 1) + (Nx + 1 - NxM); 

% Lower Legs
j = [1:NxL+1, NxM+1:Nx+1];
for i = 1:NyL
    j_ctr = (1 : Nsml) + (i-1)*Nsml;
   U(i,j) = TMPU( j_ctr ); 
end

% Middle of H
j = 1:Nx+1;
for i = NyL+1 : NyM+2
    j_ctr = NyL*Nsml + j + (i-NyL-1)*(Nx+1);
   U(i,j) = TMPU( j_ctr ); 
end

% Upper Legs
j = [1:NxL+1, NxM+1:Nx+1];
for i = NyM+3 : Ny+2
    j_ctr = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (1 : Nsml) + (i-NyM-3)*Nsml;
   U(i,j) = TMPU( j_ctr ); 
end




% --- Reshape V ---
V = zeros(Ny+1,Nx+2);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 2) + (Nx + 2 - NxM); 

% Lower Legs
j = [1:NxL+2, NxM+1:Nx+2];
for i = 1:NyL
    j_ctr = (1 : Nsml) + (i-1)*Nsml;
    V(i,j) = TMPV( j_ctr ); 
end

% Middle of H
j = 1:Nx+2;
for i = NyL+1 : NyM+1
    j_ctr = NyL*Nsml + j + (i-NyL-1)*(Nx+2);
    V(i,j) = TMPV( j_ctr ); 
end

% Upper Legs
j = [1:NxL+2, NxM+1:Nx+2];
for i = NyM+2 : Ny+1
    j_ctr = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (1 : Nsml) + (i-NyM-2)*Nsml;
    V(i,j) = TMPV( j_ctr ); 
end
