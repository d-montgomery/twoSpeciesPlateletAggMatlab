function [C,Cghst] = react_ode_solve(N,R,R0,C,dt)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Wash Channel 
i = 1:Ny;
for j = 1:NxL
    C(i+1,j+1) = C(i+1,j+1) + 1.5*dt*R(i,j) - 0.5*dt*R0(i,j);
end

% Blood Channel
for j = NxM+1:Nx
    C(i+1,j+1) = C(i+1,j+1) + 1.5*dt*R(i,j) - 0.5*dt*R0(i,j);
end

% Injury Channel
i = NyL+1:NyM;
for j = NxL+1:NxM
    C(i+1,j+1) = C(i+1,j+1) + 1.5*dt*R(i,j) - 0.5*dt*R0(i,j);
end

% % Boundary Conditions (no flux everywhere)
% C(Ny+2,:) = C(Ny+1,:); % \partial \Omega_1,3
% C(1,:) = C(2,:); % \partial \Omega_2,4
% C(NyM+1:Ny+2,NxL+2) = C(NyM+1:Ny+2,NxL+1); % \Gamma_1
% C(1:NyL+1,NxL+2) = C(1:NyL+1,NxL+1); % \Gamma_2
% C(NyM+1:Ny+2,NxM+1) = C(NyM+1:Ny+2,NxM+2); % \Gamma_3
% C(1:NyL+1,NxM+1) = C(1:NyL+1,NxM+2); % \Gamma_4
% C(NyM+1,NxL+3:NxM) = C(NyM,NxL+3:NxM); % \Gamma_5
% C(NyL,NxL+3:NxM) = C(NyL+1,NxL+3:NxM); % \Gamma_6
% C(:,1) = C(:,2); % \Gamma_7
% C(:,Nx+2) = C(:,Nx+1); % \Gamma_8
% 
% % Corners
% Cghst(1) = C(NyM+1,NxL+2);
% Cghst(2) = C(NyL+2,NxL+2);
% Cghst(3) = C(NyM+1,NxM+1);
% Cghst(4) = C(NyL+2,NxM+1);

Cghst = zeros(4,1);

% \partial \Omega_1 CN = CP
i = Ny; j = 1:NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig+1,jg) = C(ig,jg);

% \partial \Omega_2 CS = CP
i = 1; j = 1:NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig-1,jg) = C(ig,jg);

% \partial \Omega_3 CN = CP
i = Ny; j = NxM+1:Nx;
ig = i+1; jg = j+1; % Ghost nodes
C(ig+1,jg) = C(ig,jg);

% \partial \Omega_4 CS = CP
i = 1; j = NxM+1:Nx;
ig = i+1; jg = j+1; % Ghost nodes
C(ig-1,jg) = C(ig,jg);

% \Gamma_1 Boundary CE = CP
i = NyM+1:Ny; j = NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg+1) = C(ig,jg);

% \Gamma_2 Boundary CE = CP
i = 1:NyL; j = NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg+1) = C(ig,jg);

% \Gamma_3 Boundary CW = CP
i = NyM+1:Ny; j = NxM+1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg-1) = C(ig,jg);

% \Gamma_4 Boundary CW = CP
i = 1:NyL; j = NxM+1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg-1) = C(ig,jg);

% \Gamma_5 Boundary CN = CP
i = NyM; j = NxL+2:NxM-1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig+1,jg) = C(ig,jg);

% \Gamma_5 Corners CN = CP
i = NyM;
j = NxL+1;
ig = i+1; jg = j+1; % Ghost nodes
Cghst(1) = C(ig,jg);
j = NxM; jg = j+1;
Cghst(3) = C(ig,jg);

% \Gamma_6 Boundary CS = CP
i = NyL+1; j = NxL+2:NxM-1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig-1,jg) = C(ig,jg);

% \Gamma_6 Corners CS = CP
i = NyL+1; j = NxL+1;
ig = i+1; jg = j+1; % Ghost nodes
Cghst(2) = C(ig,jg);
j = NxM; jg = j+1;
Cghst(4) = C(ig,jg);

% \Gamma_7 Boundary CW = CP
i = 1:Ny; j = 1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg-1) = C(ig,jg);

% \Gamma_8 Boundary CE = CP
i = 1:Ny; j = Nx;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg+1) = C(ig,jg);

