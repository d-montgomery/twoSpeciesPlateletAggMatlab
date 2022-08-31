function h = inletBC_Pm(Plchar,WL,Wmid,Xchar,XC,N,t,Tchar)
% Unpack Various N values for H-domain
Ny = N(6);

x = XC(Ny+2,:);

% All quantities are dimensionless (except characteristic scales)
h3 = @(x) heaviside(Xchar*x - Xchar*(WL + Wmid))*Plchar; % Right Only
% h3 = @(x) Plchar - heaviside(Xchar*x - Xchar*WL)*Plchar; % Left Only

t_s = 0.001 / Tchar;
r = @(t) t/t_s + heaviside(t-t_s).*(1 - t/t_s);

h = r(t).*h3(x)/Plchar;

% Both Channels
% h = r(t) + 0*x;
