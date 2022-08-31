function [] = makeStreamlinesMobileP(N,D,U,V,XU,YU,XV,YV,XP,YP,t,B,maxB)
% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

WL = D(3); Wmid = D(4); WR = D(5); 
Hlow  = D(6); Hmid = D(7); Hup = D(8);





Ns = 75;

[U,YU,V,XV] = adjustGrid(N,U,YU,V,XV);
xu = XU(1,:); yv = YV(:,1);

[~,~,~,~,xp,yp,~] = NSGridH(Ns,Ns,WL,Wmid,WR,Hlow,Hmid,Hup,0,0);
[XPs,YPs] = meshgrid(xp(2:Ns+1),yp(2:Ns+1));

U = interp2(XU, YU, U, XPs,YPs);
V = interp2(XV, YV, V, XPs,YPs);
pcolor(XP,YP,B)
shading flat;
colormap autumn;
% set(gca,'ColorScale','log')
slines = streamline(XPs,YPs,U,V,XPs(1,:), YPs(end,1) + 0 *XPs(1,:));
set(slines,'Color','k')
axis([0 xu(end) 0 yv(end)])
colorbar
% caxis([0 maxB])

gry =  [0.7 0.7 0.7];
pos = [WL, 0, Wmid, Hlow];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')
pos = [WL, Hlow+Hmid, Wmid, Hlow+Hmid+Hup];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')

xlabel('x', 'FontSize',14)
ylabel('y', 'FontSize',14)
xlabs = [xu(1), xu(NxL+1), xu(NxM+1), xu(Nx+1)];
xticks(xlabs)

ylabs = [yv(1), yv(NyL+1), yv(NyM+1), yv(Ny+1)];
yticks(ylabs)
title(['$\mathbf{\psi(x,y,t)}$ \textbf{at} $\mathbf{t =}$',...
    num2str(t)],'Interpreter', 'latex', 'FontSize',20)

