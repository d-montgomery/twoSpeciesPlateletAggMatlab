function [C,injGhst] = ghostNodesC(N,C,BC)
% Calculates the ghost nodes in advection-diffusion equation
% Inputs:
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% C = concentration 
% xc and yc is where c(x,y,t) is stored on grid
% BC = cell array with vectors containing boundary conditions
% t = t^{n+1}

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

zy = zeros(Ny,1); zx = zeros(1,Nx+2);
C = [     zx    ;
      zy,  C , zy;
          zx     ]; % Augment for ghost nodes

injGhst = zeros(4,1);

% \partial \Omega_1 CN = -CP + 2*k1(XP,t);
i = Ny; j = 1:NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig+1,jg) = -C(ig,jg) + 2*BC(jg);

% \partial \Omega_2 CS = CP - (YP - YS) * k2(XP,t);
i = 1; j = 1:NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig-1,jg) = C(ig,jg);

% \partial \Omega_3 CN = -CP + 2*k3(XP,t);
i = Ny; j = NxM+1:Nx;
ig = i+1; jg = j+1; % Ghost nodes
C(ig+1,jg) = - C(ig,jg) + 2 * BC(jg);

% \partial \Omega_4 CS = CP - (YP - YS) * k4(XP,t);
i = 1; j = NxM+1:Nx;
ig = i+1; jg = j+1; % Ghost nodes
C(ig-1,jg) = C(ig,jg);

% \Gamma_1 Boundary CE = CP + (XE - XP) * wc_1(YP,t);
i = NyM+1:Ny; j = NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg+1) = C(ig,jg);

% \Gamma_2 Boundary CE = CP + (XE - XP) * wc_2(YP,t);
i = 1:NyL; j = NxL;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg+1) = C(ig,jg);

% \Gamma_3 Boundary CW = CP - (XP - XW) * wc_3(YP,t);
i = NyM+1:Ny; j = NxM+1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg-1) = C(ig,jg);

% \Gamma_4 Boundary CW = CP - (XP - XW) * wc_4(YP,t);
i = 1:NyL; j = NxM+1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg-1) = C(ig,jg);

% \Gamma_5 Boundary CN = CP + (YN - YP)*wc_5(XP,t);
i = NyM; j = NxL+2:NxM-1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig+1,jg) = C(ig,jg);

% \Gamma_5 Corners CN = CP + (YN - YP)*wc_5(XP,t); 
i = NyM;
j = NxL+1;
ig = i+1; jg = j+1; % Ghost nodes
injGhst(1) = C(ig,jg);
j = NxM; jg = j+1;
injGhst(3) = C(ig,jg);

% \Gamma_6 Boundary CS = CP - (YP - YS) * wc_6(XP,t);
i = NyL+1; j = NxL+2:NxM-1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig-1,jg) = C(ig,jg);

% \Gamma_6 Corners CS = CP - (YP - YS) * wc_6(XP,t);
i = NyL+1; j = NxL+1;
ig = i+1; jg = j+1; % Ghost nodes
injGhst(2) = C(ig,jg);
j = NxM; jg = j+1;
injGhst(4) = C(ig,jg);

% \Gamma_7 Boundary CW = CP - (XP - XW) * wc_7(YP,t);
i = 1:Ny; j = 1;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg-1) = C(ig,jg);

% \Gamma_8 Boundary CE = CP + (XE - XP) * wc_8(YP,t);
i = 1:Ny; j = Nx;
ig = i+1; jg = j+1; % Ghost nodes
C(ig,jg+1) = C(ig,jg);