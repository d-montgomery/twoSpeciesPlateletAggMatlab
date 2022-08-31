function k = k_adh_fxn(XC,YC,Xchar)
%XC, YC should have dimensions cm

% act_rate = (2e-10/6.022); % standard slow rate
% act_rate = (2e-8/6.022); % faster rate for test
act_rate = 4;

% Center of elipse on bottom
xm = 3.5e-2 / 2 / Xchar;
ym  = .75e-2 / Xchar; % Bottom of channel

% Radius
a = 8.5e-4 / Xchar; % major axis (width = 2a)
b = 8.5e-4 / Xchar; % minor axis (height = 2b) or 8.5
% a = 0.7e-2 / Xchar; % major axis (width = 2a)
% b = 10e-4 / Xchar; % minor axis (height = 2b) or 8.5

% xi(x,y) is elipse function , smoothes it a bit and scales it by act_rate
xi = sqrt(((XC - xm)./a).^2 + ((YC - ym)./b).^2) - 1;
k = 0.5*(1 - erf(1e2*xi))*act_rate; % cm^3 / s


% % Center of elipse on Top
% xm = 3.5e-2 / 2 / Xchar;
% ym  = .95e-2 / Xchar; % Top of channel
% 
% % Radius
% a = 0.70e-2 / Xchar; % major axis (width = 2a)
% b = 10e-4 / Xchar; % minor axis (height = 2b)
% 
% % xi(x,y) is elipse function , smoothes it a bit and scales it by act_rate 
% xi = sqrt(((XC - xm)./a).^2 + ((YC - ym)./b).^2) - 1;
% 
% k = k + 0.5*(1 - erf(1e2*xi))*act_rate; % cm^3 / s
% 
% 
