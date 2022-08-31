function U0 = initial_U(N,QM,Hlow,Hmid,Xchar,Uchar,XU,YU)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3); Ny = N(6);

xi = 6*QM/(Xchar*Hmid)^3;
Iy = @(y) min( xi*(Xchar * y - Xchar*Hlow - Xchar*Hmid/2).^2 ...
          - xi*(Xchar*Hmid/2)^2, 0)/Uchar ; % Top Right
ICu = @(x,y) Iy(y);
U0 = ICu(XU,YU);
U0(1:Ny+1,1:NxL) = 0;
U0(1:Ny+1,NxM+1:Nx+1) = 0;