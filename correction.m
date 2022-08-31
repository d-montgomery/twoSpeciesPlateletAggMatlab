function[U,V,P] = correction(U0,V0,P,PHI,N,dt,xp,yp,BC_in)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

U = 0*U0; 
V = 0*V0;

% update Pressure
P = P + PHI;
% Zero out interior of H
P(1:NyL,NxL+1:NxM) = 0;
P(NyM+1:Ny,NxL+1:NxM) = 0;

% update U - Wash Channel
i = 2 : Ny+1;
j = 2 : NxL;
U(i,j) = U0(i,j) - dt * ( PHI(i-1,j) - PHI(i-1,j-1) )./( xp(j+1)-xp(j) );
U(1,1:NxL+1) = U(2,1:NxL+1); % du/dy = 0 on bottom boundary
U(Ny+2,j) = 0; % u = 0 on top boundary

% update U - Injury Channel
i = NyL+2 : NyM+1; %%% NyM
j = NxL+1 : NxM+1;
U(i,j) = U0(i,j) - dt * ( PHI(i-1,j) - PHI(i-1,j-1) )./( xp(j+1)-xp(j) );

% update U - Blood Channel
i = 2 : Ny+1;
j = NxM+2 : Nx;
U(i,j) = U0(i,j) - dt * ( PHI(i-1,j) - PHI(i-1,j-1) )./( xp(j+1)-xp(j) );
U(1,NxM+1:Nx+1) = U(2,NxM+1:Nx+1); % du/dy = 0 on bottom boundary
U(Ny+2,j) = 0; % u = 0 on top boundary

% update V - Wash Channel
i = 2 : Ny;
j = 2 : NxL+1;
V(i,j) = V0(i,j) - dt * ( PHI(i,j-1) - PHI(i-1,j-1) )./( yp(i+1)-yp(i) )';
V(1,1:NxL+2) = V(2,1:NxL+2); % dv/dy = 0 on bottom boundary
V(Ny+1,j) = BC_in(j); % v = inlet condition on top

% update V - Injury Channel
i = NyL+2 : NyM;
j = NxL+2 : NxM+1;
V(i,j) = V0(i,j) - dt * ( PHI(i,j-1) - PHI(i-1,j-1) )./( yp(i+1)-yp(i) )';


% update V - Blood Channel
i = 2 : Ny;
j = NxM+2 : Nx+1;
V(i,j) = V0(i,j) - dt * ( PHI(i,j-1) - PHI(i-1,j-1) )./( yp(i+1)-yp(i) )';
V(1,NxM+1:Nx+2) = V(2,NxM+1:Nx+2); % dv/dy = 0 on bottom boundary
V(Ny+1,j) = BC_in(j); % v = inlet condition on top