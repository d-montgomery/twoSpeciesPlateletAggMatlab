function [Pb, PbInjGhst] = Pb_soln(N,XC,YC,Pb_max,t,Plchar,Xchar)
%XC, YC should have dimensions cm
%Pb_max should have dimensions platelets/vol

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); 
NyL = N(4); NyM = N(5); 

% Pb = f(t)g(x,y)
% Midpoints of circle
xm = 3.5e-2 / 2;
ym  = .75e-2; % Bottom of channel
% ym = .95e-2; % Top of channel

% Radius
r = 8.5e-4;

% xi(x,y) circular function 
xi = sqrt((XC - xm).^2 + (YC - ym).^2) - r;

% g(x,y) function  
g = 0.5*(1 - erf(1e4*xi)); % smooth edges

% f(t) function
f = @(t) Pb_max /Plchar;

Pb = f(t)*g;

PbInjGhst = [Pb(NyM+2,NxL+2); Pb(NyL+1,NxL+2); ...
             Pb(NyM+2,NxM+1); Pb(NyL+1,NxM+1)];
