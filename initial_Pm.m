function [Pm, PmInjGhst] = initial_Pm(N,XC,YC)
% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2);
NyL = N(4); NyM = N(5); 

Pm = 0*XC.*YC;
PmInjGhst = [Pm(NyM+2,NxL+2); Pm(NyL+1,NxL+2); ...
             Pm(NyM+2,NxM+1); Pm(NyL+1,NxM+1)];