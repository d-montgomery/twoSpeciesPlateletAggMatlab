function [] = makeVectorfield_PbPm(D,N,U,V,XU,YU,XV,YV,XP,YP,Pb,Pm,t,Pmax)
% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Rescale to micro meters
D = 1e4*D;
XU = 1e4*XU; YU = 1e4*YU;
XV = 1e4*XV; YV = 1e4*YV;
XP = 1e4*XP; YP = 1e4*YP;
Pb = 1e-3*Pb; Pm = 1e-3*Pm;
Pmax = 1e-3*Pmax;

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

% Create pcolor of mobile platelets 
ax1 = axes;
pm = pcolor(XP,YP,Pm); shading flat;
hold on

% Create pcolor of bound platelets 
Pb(abs(Pb)<1e-8)=0;
ax2 = axes;
if max(max(Pb)) > 1e-8 % Don't plot Pb until Pb > 1e-8
    pb = pcolor(XP,YP,Pb); shading flat;
    alpha(pb,'color')
    alpha(pb, 'scaled')
    hold on
end


% Create vectorfield plot
ax3 = axes;
veclines = quiver(XPs, YPs, U,V);
set(veclines,'Color','k')
axis([XP(1,1) XP(1,end) YP(1,1) YP(end,1)])

axis([XP(1,1) XP(1,end) YP(1,1) YP(end,1)])
gry =  [0.7 0.7 0.7];
pos = [WL, YP(1,1), Wmid, Hlow];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')
pos = [WL, Hlow+Hmid, Wmid, YP(end,end)];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')
% % Boarder around figure
% pos = [XP(1,1), YP(1,1), XP(end,end), YP(end,end)];
% rectangle('Position', pos,'EdgeColor','k','LineWidth',1)


% Give each Platelet its own colormap
% position = [left, bottom, width, height]
colormap(my_colormaps('cmp_bluetored'))
% colormap(ax1,'autumn')
colormap(ax2,'spring')


% Align the axes
% set(gcf,'Position',[100 100 500 500])
% set([ax1,ax2,ax3],'Position',[.17 .11 .685 .815]);
set([ax1,ax2,ax3],'Position',[.18 .11 .65 .815]);
cb1 = colorbar(ax1,'Position',[.07 .11 .03 .815]);
cb2 = colorbar(ax2,'Position',[.86 .11 .03 .815]);
ylabel(cb1,'Mobile Platelets per \muL','FontSize',16)
% Can't make right colorbar label until end?????
% ylabel(cb2,'Bound Platelets per \muL','FontSize',16,'Rotation',270)
% cb2.Label.Position(1) = 3.5;
% cb2.Label.Position = [pos(1)+3, pos(2)+390];
caxis(ax1,[0,max(max(Pm))])
caxis(ax2,[0,max(1e4,max(max(Pb)))])
% ax1.Visible = 'off';
ax2.Visible = 'off';
ax3.Visible = 'off';
xlabs = [XP(1,1), XU(1,NxL+1), XU(1,NxM+1), XP(1,Nx)];
ylabs = [YP(1,1), YV(NyL+1,1), YV(NyM+1,1), YP(Ny,1)];
ax1.XTick = xlabs;
ax1.YTick = ylabs;
xlabs = [0, XU(1,NxL+1), XU(1,NxM+1), XU(1,Nx+1)];
ylabs = [0, YV(NyL+1,1), YV(NyM+1,1), YV(Ny+1,1)];
ax1.XLabel.String = 'x (\mum)';
ax1.YLabel.String = 'y (\mum)';
ax1.XTickLabel = strsplit(num2str(xlabs));
ax1.YTickLabel = strsplit(num2str(ylabs));
ax1.Title.String = ['Mobile and Bound Platelets at t = ',num2str(t),' s'];
ax1.Title.FontSize = 16;

% Label for right colorbar
ylabel(cb2,'Bound Platelets per \muL','FontSize',16,'Rotation',270)
cb2.Label.Position(1) = 4;
set(cb2.Label,'FontSize',16)






