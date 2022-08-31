function [] = makeVectorfieldBrinkman(D,N,U,V,XU,YU,XV,YV,XP,YP,B,maxB,t)
% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Get coarse grid for plotting vector field 
WL = D(3); Wmid = D(4); WR = D(5); 
Hlow  = D(6); Hmid = D(7); Hup = D(8);

[U,YU,V,XV] = adjustGrid(N,U,YU,V,XV);
xu = XU(1,:); yv = YV(:,1);

Nxs = 50;
Nys = Nxs;
[~,~,~,~,xps,yps,~] = NSGridH(Nxs,Nys,WL,Wmid,WR,Hlow,Hmid,Hup,0,0);
[XPs,YPs] = meshgrid(xps(2:Nxs+1),yps(2:Nys+1)); %Ignore ghost cells

% Interpolate U and V to cell centers
U = interp2(XU,YU,U,XPs,YPs);
V = interp2(XV,YV,V,XPs,YPs);

% Create pcolor of brinkman term
pcolor(XP,YP,B); shading flat;
colormap autumn;
hold on

% Create vectorfield plot
veclines = quiver(XPs, YPs, U,V);
set(veclines,'Color','k')
axis([0 xu(end) 0 yv(end)])
colorbar
caxis([min(min(B)) max(max(B))])

axis([0 XU(1,end) 0 YV(end,1)])
gry =  [0.7 0.7 0.7];
pos = [WL, 0, Wmid, Hlow];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')
pos = [WL, Hlow+Hmid, Wmid, Hlow+Hmid+Hup];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')

xlabel('x', 'FontSize',14)
ylabel('y', 'FontSize',14)
xlabs = [XU(1,1), XU(1,NxL+1), XU(1,NxM+1), XU(1,Nx+1)];
xticks(xlabs)

ylabs = [YV(1,1), YV(NyL+1,1), YV(NyM+1,1), YV(Ny+1,1)];
yticks(ylabs)
title(['Mobile Platelets at t = ',num2str(t)],'FontSize',16);
