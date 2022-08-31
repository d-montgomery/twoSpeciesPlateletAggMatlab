function a = alpha_fxn(T,Bchar)
% Frictional Resistance (alpha function)
C_CK = 1e8; % 1 / cm^2
a = C_CK/Bchar * (0.6*T).^2 ./ (1-0.6*T).^3;