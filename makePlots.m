function [] = makePlots(U,V,P,XU,YU,XV,YV,XP,YP,t,D,N)
WL = D(3); Wmid = D(4); WR = D(5); 
Hlow  = D(6); Hmid = D(7); Hup = D(8);

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);


gry =  [0.7 0.7 0.7];
posL = [WL,0, Wmid, Hlow];
posU = [WL,Hlow+Hmid ,Wmid,Hup];

axisvec = [0, WL + Wmid + WR, 0, Hlow + Hmid + Hup];

labelFont = 16;
        
% Plot the U solution
subplot(1,3,1)
% surf(XU,YU,U)
pcolor(XU,YU,U)
colorbar()
shading flat
rectangle('Position', posL, 'FaceColor',gry,'EdgeColor',gry)
rectangle('Position', posU, 'FaceColor',gry,'EdgeColor',gry)
xlabel('x', 'FontSize', labelFont)
ylabel('y', 'FontSize', labelFont)
zlabel('U', 'FontSize', labelFont)
title('U', 'FontSize', labelFont)
xlabs = [0, WL, WL+Wmid, WL + Wmid + WR];
xticks(xlabs)
ylabs = [0, Hlow, Hlow + Hmid, Hlow + Hmid + Hup];
yticks(ylabs)
axis(axisvec)

% Plot the V Solution
subplot(1,3,2)
% surf(XV,YV,V))
pcolor(XV,YV,V)
colorbar()
shading flat
rectangle('Position', posL, 'FaceColor',gry,'EdgeColor',gry)
rectangle('Position', posU, 'FaceColor',gry,'EdgeColor',gry)
xlabel('x', 'FontSize', labelFont)
ylabel('y', 'FontSize', labelFont)
zlabel('V', 'FontSize', labelFont)
title('V', 'FontSize', labelFont)
xlabs = [0, WL, WL+Wmid, WL + Wmid + WR];
xticks(xlabs)
ylabs = [0, Hlow, Hlow + Hmid, Hlow + Hmid + Hup];
yticks(ylabs)
axis(axisvec)

% Plot the P Solution
subplot(1,3,3)
pcolor(XP,YP,P)
mxm1 = max(max(P(:,1:NxL))); mnm1 = min(min(P(:,1:NxL)));
mxm2 = max(max(P(:,NxM+1:Nx))); mnm2 = min(min(P(:,NxM+1:Nx)));
colorbar()
caxis([min(mnm1,mnm2), max(mxm1,mxm2)])

shading flat
rectangle('Position', posL, 'FaceColor',gry,'EdgeColor',gry)
rectangle('Position', posU, 'FaceColor',gry,'EdgeColor',gry)
xlabel('x', 'FontSize', labelFont)
ylabel('y', 'FontSize', labelFont)
zlabel('P', 'FontSize', labelFont)
title('P', 'FontSize', labelFont)
xlabs = [0, WL, WL+Wmid, WL + Wmid + WR];
xticks(xlabs)
ylabs = [0, Hlow, Hlow + Hmid, Hlow + Hmid + Hup];
yticks(ylabs)
axis(axisvec)

sgtitle(['$\vec u$ and $p$ at $t = $ ', num2str(t)],...
         'Interpreter', 'latex', 'FontSize', 24) 

