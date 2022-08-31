function [Bp,dx,dy] = Bp_Ucoeff(i,j,XU,XV,YV,B)
    dy = YV(i,j)-YV(i-1,j);
    dx = XV(i,j+1) - XV(i,j);
    XP = XU(i,j);
    xe = XV(i,j+1); 
    xw = XV(i,j);
    Be = B(i-1,j); Bw = B(i-1,j-1);
%     Bp = (Be*(xe-XP) + Bw*(XP-xw))/dx;
    Bp = ( Be*(XP - xw)-Bw*(XP - xe) ) / dx; 
end