function [] = makePcolorsPlatelets(D,XP,YP,t,Pm,maxPm,Pb,maxPb)
% Unpack Various N values for H-domain

WL = D(3); Wmid = D(4); WR = D(5); 
Hlow  = D(6); Hmid = D(7); Hup = D(8);

% Mobile Species Plot
subplot(1,2,1)
pcolor(XP,YP,Pm)
shading flat;
colormap autumn;
axis([0 WL+Wmid+WR 0 Hlow+Hmid+Hup])
colorbar
% caxis([0 maxPm])

gry =  [0.7 0.7 0.7];
pos = [WL, 0, Wmid, Hlow];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')
pos = [WL, Hlow+Hmid, Wmid, Hlow+Hmid+Hup];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')

xlabel('x', 'FontSize',14)
ylabel('y', 'FontSize',14)
xlabs = [0, WL, WL+Wmid, WL+Wmid+WR];
xticks(xlabs)

ylabs = [0, Hlow, Hlow+Hmid, Hlow+Hmid+Hup];
yticks(ylabs)
title(['Mobile Platelets at t =',...
    num2str(t)], 'FontSize',20)


% Bound Species Plot
subplot(1,2,2)
pcolor(XP,YP,Pb)
shading flat;
colormap autumn;
axis([0 WL+Wmid+WR 0 Hlow+Hmid+Hup])
colorbar
% caxis([0 maxPb])

gry =  [0.7 0.7 0.7];
pos = [WL, 0, Wmid, Hlow];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')
pos = [WL, Hlow+Hmid, Wmid, Hlow+Hmid+Hup];
rectangle('Position', pos, 'FaceColor',gry,'EdgeColor','k')

xlabel('x', 'FontSize',14)
ylabel('y', 'FontSize',14)
xlabs = [0, WL, WL+Wmid, WL+Wmid+WR];
xticks(xlabs)

ylabs = [0, Hlow, Hlow+Hmid, Hlow+Hmid+Hup];
yticks(ylabs)
title(['Bound Platelets at t =',...
    num2str(t)], 'FontSize',20)

