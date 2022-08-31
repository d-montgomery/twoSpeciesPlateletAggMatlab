function [Bp,dx,dy] = Bp_Vcoeff(i,j,XU,YU,YV,B)
dy = YU(i+1,j)-YU(i,j);
YP = YV(i,j);
dx = XU(i,j) - XU(i,j-1);
yn = YU(i+1,j);
ys = YU(i,j);
Bn = B(i,j-1);
Bs = B(i-1,j-1);
Bp = (Bn*(YP-ys) - Bs*(YP-yn))/dy;
end