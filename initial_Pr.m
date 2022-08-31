function [P0] = initial_Pr(N,D,P_io,Xchar,Pchar,XP,YP)

% unpack dimensions of H
WL = D(3); Wmid = D(4); WR = D(5); 
Hlow  = D(6); Hmid = D(7); Hup = D(8);

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Initial Conditions 
Pw_in = P_io(1); Pw_o = P_io(2);
Pb_in = P_io(3); Pb_o = P_io(4);

p_w = @(y) ( Pchar*(Pw_in - Pw_o)/(Xchar*(Hlow + Hmid +Hup))...
             * Xchar*y ) + Pchar*Pw_o;
p_b = @(y) ( Pchar*(Pb_in - Pb_o)/(Xchar*(Hlow + Hmid +Hup))...
             * Xchar*y ) + Pchar*Pb_o;

A = Xchar*[1/Xchar, WL,      Hlow; 
           1/Xchar, WL,      Hlow+Hmid; 
           1/Xchar, WL+Wmid, Hlow;
           1/Xchar, WL+Wmid, Hlow+Hmid];
b = [p_w(Hlow);
     p_w(Hlow+Hmid);
     p_b(Hlow);
     p_b(Hlow+Hmid)];
 
C = (A'*A) \ (A'*b);

p_in = @(x,y) (C(1) + C(2)*Xchar*x + C(3)*Xchar*y);

Pbw = @(x,y) heaviside(Xchar*WL-Xchar*x).*p_w(y)... 
               + heaviside(Xchar*x-Xchar*(WL+Wmid)).*p_b(y);
P0 = Pbw(XP,YP);
P0(NyL+1:NyM,NxL+1:NxM) = p_in(XP(NyL+1:NyM,NxL+1:NxM),YP(NyL+1:NyM,NxL+1:NxM));
P0 = P0/Pchar;