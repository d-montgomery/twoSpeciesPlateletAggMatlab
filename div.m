function D = div(U,V,XU,YV, Nx, Ny)

% D = zeros(Ny,Nx);
% for i = 1:Ny
%     for j = 1:Nx
%        dx = xu(j+1) - xu(j);
%        dy = yv(i+1) - yv(i);
%        
%        ue = U(i+1,j+1); uw = U(i+1,j);
%        vn = V(i+1,j+1); vs = V(i,j+1);
%        D(i,j) = (ue - uw)./dx + (vn - vs)./dy;
%     end
% end

dx = XU(2:Ny+1,2:Nx+1) - XU(2:Ny+1,1:Nx);
dy = YV(2:Ny+1,2:Nx+1) - YV(1:Ny,2:Nx+1);

ue = U(2:Ny+1 , 2:Nx+1); uw = U(2:Ny+1 , 1:Nx);
vn = V(2:Ny+1 , 2:Nx+1); vs = V(1:Ny , 2:Nx+1);

D = (ue - uw)./dx + (vn - vs)./dy;
