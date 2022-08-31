function V0 = initial_V(Q1,Q3,WL,Wmid,WR,Xchar,Uchar,XV)
alpha = 6*Q1/(Xchar*WL)^3;
beta = 6*Q3/(Xchar*WR)^3;
W = WL + Wmid + WR;
h1 = @(x) min( alpha*(Xchar * x-Xchar*WL/2).^2 ...
                 - alpha*(Xchar*WL/2)^2 , 0)/Uchar; % Top Left
h3 = @(x) min( beta*(Xchar * x - Xchar*W + Xchar*WR/2).^2 ...
               - beta*(Xchar*WR/2)^2, 0)/Uchar ;
V0 = h1(XV) + h3(XV);