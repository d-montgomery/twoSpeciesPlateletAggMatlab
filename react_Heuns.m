function [Pm,Pb] = react_Heuns(N,Pm,Pb,k_adh,dt,Pmax,Plchar,Rchar)
% React_Pb(Pm,Pb,k_adh,Pmax)
% React_Pb(Pm,Pb,k_adh,Pmax)
% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

Rm = @(Pm,Pb,k_adh) React_Pb(Plchar*Pm,Plchar*Pb,k_adh,Pmax,Rchar);
Rb = @(Pm,Pb,k_adh) React_Pm(Plchar*Pm,Plchar*Pb,k_adh,Pmax,Rchar);

% Wash Channel 
i = 1:Ny;
for j = 1:NxL
    [Pm(i,j),Pb(i+1,j+1)] = heuns(Rm,Rb,Pm(i,j),Pb(i+1,j+1),...
                                      k_adh(i,j),dt);
end

% Blood Channel
for j = NxM+1:Nx
    [Pm(i,j),Pb(i+1,j+1)] = heuns(Rm,Rb,Pm(i,j),Pb(i+1,j+1),...
                                      k_adh(i,j),dt);
end

% Injury Channel
i = NyL+1:NyM;
for j = NxL+1:NxM
    [Pm(i,j),Pb(i+1,j+1)] = heuns(Rm,Rb,Pm(i,j),Pb(i+1,j+1),...
                                      k_adh(i,j),dt);
end

end

function [Pm,Pb] = heuns(Rm,Rb,Pm,Pb,k_adh,dt)
    k1m = Rm(Pm,Pb,k_adh);
    k1b = Rb(Pm,Pb,k_adh);
    k2m = Rm(Pm + dt*k1m, Pb + dt*k1b, k_adh);
    k2b = Rb(Pm + dt*k1m, Pb + dt*k1b, k_adh);
    Pm = Pm + 0.5*dt*(k1m+k2m);
    Pb = Pb + 0.5*dt*(k1b+k2b);
end
