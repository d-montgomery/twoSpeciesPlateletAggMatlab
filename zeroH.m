function [U,V,P] = zeroH(U,V,P,NxL,NyL,Nxm,Nym,Ny)

% Force zeros in middle to form H (doesn't zero ghost nodes)
U(1:NyL,NxL+2:Nxm) = 0;
U(Nym+3:Ny+2,NxL+2:Nxm) = 0;
V(1:NyL,NxL+3:Nxm) = 0;
V(Nym+2:Ny+1,NxL+3:Nxm) = 0;

% Zeros out ghost nodes
P(1:NyL,NxL+1:Nxm) = 0;
P(Nym+1:Ny,NxL+1:Nxm) = 0;

end