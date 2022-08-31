function [XU,YU,XV,YV,xp,yp,N,mdx,mdy] = ...
             NSGridH(Nx,Ny,Wleft,Wmid,Wright,Hlow,Hmid,Hup,GRID,PLOT)

% Inputs: Nx and Ny are number of cells in x and y direction
%         Lx and Ly are the width and height of H
%         Ix and Iy are the width and height of injury channel
%         GRID is a flag for different types of grids
%         PLOT = 1 produces plot of staggered grid

% Determine overall width and height of Domain
W = Wleft + Wmid + Wright;
H = Hlow + Hmid + Hup;

% Determine number of grid points in each region
if GRID == 0
    % Number of x-grid points for each region
    NxL = round(Nx*Wleft/W); % x pts in wash channel
    NxR = round(Nx*Wright/W); % x pts in blood channels
    NxMid = round(Nx*Wmid/W); % Nx in Injury Channel
    Nx = NxMid + NxL + NxR; % Total # of points in x-direction

    % Number of y-grid points for each region
    NyL = round(Hlow/H*Ny); % Lower Legs of H
    NyMid = round(Hmid/H*Ny); % Injury Channel
    NyUp = round(Hup/H*Ny); % Upper Legs of H
    Ny = NyL + NyMid + NyUp; % Total # of points in y-direction
    NxM = NxL + NxMid;
    NyM = NyL + NyMid;
    N = [NxL, NxM, Nx, NyL, NyM, Ny];
else
    % percentage of x-pts and y-pts in injury channel
    Sx = 0.55; 
    Sy = 0.3;

    
    % Number of x-grid points for each region
    NxMid = round(Sx*Nx);
    WashPercent = Wleft/(W - Wmid);
    NxL = round((1-Sx)*Nx*WashPercent);
    BloodPercent = Wright/(W - Wmid);
    NxR = round((1-Sx)*Nx*BloodPercent);
    
    % percentage of y-pts in injury channel
    NyMid = round(Sy*Ny);
    NyL = round(Hlow/(H-Hmid)*(1-Sy)*Ny); % Injury Channel
    NyUp = round(Hup/(H-Hmid)*(1-Sy)*Ny); % Upper Legs of H
    
    Nx = NxL + NxMid + NxR;
    Ny = NyL + NyMid + NyUp;
    NxM = NxL + NxMid;
    NyM = NyL + NyMid;

    N = [NxL, NxM, Nx, NyL, NyM, Ny];
end

if GRID == 0 % "Uniform" Grid
    
    % Wash Channel - x
    dx = Wleft / NxL;
    xW = 0 : dx : Wleft;
    
     % Injury Channel - x
    dx = Wmid / NxMid;
    xI = Wleft : dx : Wleft + Wmid;
    
    % Blood Channel - x
    dx = Wright / NxR;
    xB = Wleft + Wmid : dx : W;
 
    
    % Lower - y
    dy = Hlow / NyL;
    yLow = 0 : dy : Hlow; 
 
    
    % Injury - y
    dy = Hmid / NyMid;
    yI = Hlow : dy : Hlow + Hmid;
    
    % Upper - y
    dy = Hup / NyUp;
    yUp = Hlow + Hmid : dy : H;
    
    % Build x and y vectors for grid
    x = [xW, xI(2:NxMid+1), xB(2:NxR+1)];
    y = [yLow, yI(2:NyMid+1), yUp(2:NyUp+1)];
    
else 
    
    % Injury Channel - x
    if GRID == 1
    dx = Wmid / NxMid;
    xI = Wleft : dx : Wleft + Wmid;
    else
        dxh = 0.75*(Wmid / NxMid);
        dx = (Wmid - 2*dxh) / (NxMid-2);
        xI = [Wleft, Wleft+dxh:dx:(Wleft+Wmid-dxh), Wleft + Wmid];     
    end
    
    % Blood Channel - x
    i = 2:NxL+1;
    p = log(Wleft/(xI(NxMid+1)-xI(NxMid)))/log(NxL) - 1;
    dxb = (i-1).^p * Wleft / NxL^(p+1);
    xB(1) = Wleft+Wmid;
    xB(i) = Wleft + Wmid + (i-1).*dxb;

    % Wash Channel - x
    dxw = fliplr(xB(2:NxL+1) - xB(1:NxL));
    xW = xB;
    xW(1) = 0;
    for i = 2:NxL+1
        xW(i) = xW(i-1) + dxw(i-1);
    end

    % Injury - y
    if GRID == 1
        dy = Hmid / NyMid;
        yI = Hlow : dy : Hlow + Hmid;
    else 
        dyh = 0.75 * Hmid / NyMid;
        dy = (Hmid - 2*dyh) / (NyMid-2);
        yI = [Hlow, (Hlow+dyh):dy:(Hlow+Hmid-dyh), Hlow+Hmid]; 
    end
    
    % Upper - y
    i = 2:NyUp+1;
    p = log(Hup/(yI(NyMid+1)-yI(NyMid)))/log(NyUp) - 1;
    dyu = (i-1).^p * Hup / NyUp^(p+1);
    yU(1) = Hlow+Hmid;
    yU(i) = Hlow + Hmid + (i-1).*dyu;
    
    % Lower - y
    if Hlow > 0 
        i = 2:NyL+1;
        p = log(Hlow/(yI(NyMid+1)-yI(NyMid)))/log(NyL) - 1;
        dyu = (i-1).^p * Hlow / NyL^(p+1);
        yLtmp(1) = 0;
        yLtmp(i) = 0 + (i-1).*dyu;

        dyl = fliplr(yLtmp(2:NyL+1) - yLtmp(1:NyL));
        yL(NyL+1) = Hlow;
        for i = NyL:-1:1
            yL(i) = yL(i+1) - dyl(i);
        end
    else
        yL = 0;
    end
    
    x = [xW, xI(2:NxMid), xB];
    y = [yL, yI(2:NyMid), yU];
end

% Length of each cell (for nonuniform grid)
dx = x(2:Nx+1) - x(1:Nx);
dy = y(2:Ny+1) - y(1:Ny);
mdx = min(dx);
mdy = min(dy);

% --- Location of All Points in a Rectangle ---
% Location of pressure (mid points of cells)
xp = [x(1)-dx(1)/2, x(1:Nx) + dx/2, x(Nx+1)+dx(Nx)/2];
yp = [y(1)-dy(1)/2, y(1:Ny) + dy/2, y(Ny+1)+dy(Ny)/2];

% Location of u for rectangle
xu = x;
yu = [y(1) , y(1:Ny) + dy/2, y(Ny+1) ];

% Location of v for rectangle
xv = [x(1), x(1:Nx) + dx/2, x(Nx+1)];
yv = y;

% Make Meshgrid of points for u and v
[XU,YU] = meshgrid(xu,yu);
[XV,YV] = meshgrid(xv,yv);
[XP,YP] = meshgrid(xp,yp);

% Adjust for boundaries of H
YU(NyL+1, NxL+2:NxM) = YV(NyL+1,NxL+2:NxM);
YU(NyM+2, NxL+2:NxM) = YV(NyM+1,NxL+2:NxM);
XV([1:NyL, NyM+2:Ny+1], NxL+2) = XU([1:NyL, NyM+2:Ny+1],NxL+1);
XV([1:NyL, NyM+2:Ny+1], NxM+1) = XU([1:NyL, NyM+2:Ny+1],NxM+1);



if PLOT == 1

    % Marker size
    M = 15;

    fig1 = figure;
    pos_fig1 = [0 0 750 600];
    set(fig1,'Position',pos_fig1)
    hold on
    axis([x(1) - 1.75*dx(1), x(Nx+1)+ 1.75*dx(Nx), ...
          y(1) - 1.75*dy(1), y(Ny+1)+ 1.75*dy(Ny)])
    
    % Plot vertical and horizontal grid lines
    for i = 1: Nx+1
        plot([x(i), x(i)], [y(1), y(Ny+1)], 'r');
    end
    for i = 1:Ny+1
        plot([x(1), x(Nx+1)], [y(i), y(i)], 'r');
    end
    
    % plot rectangles (covers up grid lines to make H)
    rectangle('Position',[Wleft,0,Wmid,Hlow], 'FaceColor', 'w', ...
                'EdgeColor', 'w', 'LineWidth',2.0)
    rectangle('Position',[Wleft,Hmid+Hlow,Wmid,Hup], 'FaceColor', 'w', ...
                'EdgeColor', 'w', 'LineWidth',2.0)
    

    % Plot Outline of H
    lw = 3; % line width
    plot([0,0],[0,H], 'k', 'LineWidth', lw) % left vert
    plot([W,W],[0,H], 'k', 'LineWidth', lw) % right vert
    plot([0,Wleft],[0,0], 'k', 'LineWidth', lw) % bottom W
    plot([Wleft+Wmid,W],[0,0], 'k', 'LineWidth', lw) % bottom B
    plot([0,Wleft],[H,H], 'k', 'LineWidth', lw) % top W
    plot([Wleft+Wmid,W],[H,H], 'k', 'LineWidth', lw) % top B

    plot([Wleft,Wleft+Wmid],[Hlow,Hlow], 'k', 'LineWidth', lw) % bottom I
    plot([Wleft,Wleft+Wmid],[Hlow+Hmid,Hlow+Hmid], 'k', 'LineWidth', lw) % top I
    plot([Wleft,Wleft],[0,Hlow], 'k', 'LineWidth', lw) % W right lower vert
    plot([Wleft,Wleft],[Hlow+Hmid,H], 'k', 'LineWidth', lw) % W right upper vert
    plot([Wleft+Wmid,Wleft+Wmid],[0,Hlow], 'k', 'LineWidth', lw) % B left lower vert
    plot([Wleft+Wmid,Wleft+Wmid],[Hlow+Hmid,H], 'k', 'LineWidth', lw) % B left upper vert  
    
    % Plot Wash Channel Markers
    plot(XP(:,1:NxL+2), YP(:,1:NxL+2),'b.', 'MarkerSize', M, 'DisplayName', 'Pressure');
    plot(XU(:,1:NxL+1), YU(:,1:NxL+1),'bs', 'MarkerSize', M, 'DisplayName', 'u - Velocity');
    plot(XV(:,1:NxL+2), YV(:,1:NxL+2), 'b^', 'MarkerSize', M, 'DisplayName', 'v - Velocity');
    
    % Plot Injury Channel Markers
    plot(XP(NyL+1:NyM + 2, NxL+2:NxM +1), YP(NyL+1:NyM + 2, NxL+2:NxM +1),'b.', 'MarkerSize', M, 'DisplayName', 'Pressure');
    plot(XU(NyL+1:NyM+2,NxL+2:NxM), YU(NyL+1:NyM+2, NxL+2:NxM),'bs', 'MarkerSize', M, 'DisplayName', 'u - Velocity');
    plot(XV(NyL+1:NyM+1,NxL+3:NxM), YV(NyL+1:NyM+1,NxL+3:NxM), 'b^', 'MarkerSize', M, 'DisplayName', 'v - Velocity');
    
    % Plot Blood Channel Markers
    plot(XP(:,NxM +1: Nx + 2), YP(:,NxM +1: Nx + 2),'b.', 'MarkerSize', M, 'DisplayName', 'Pressure');
    plot(XU(:,NxM+1:Nx+1), YU(:,NxM+1:Nx+1),'bs', 'MarkerSize', M, 'DisplayName', 'u - Velocity');
    plot(XV(:,NxM+1:Nx+2), YV(:,NxM+1:Nx+2), 'b^', 'MarkerSize', M, 'DisplayName', 'v - Velocity');
    
% plot(XP, YP,'b.', 'MarkerSize', M, 'DisplayName', 'Pressure');
% plot(XU, YU,'bs', 'MarkerSize', M, 'DisplayName', 'u - Velocity');
% plot(XV, YV, 'b^', 'MarkerSize', M, 'DisplayName', 'v - Velocity');
end
end