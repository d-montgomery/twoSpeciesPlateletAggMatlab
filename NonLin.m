function[NLU,NLV] = NonLin(U,V,N,XU,YU,XV,YV)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Compute NLU for x-momentum equation --------------------------------

NLU = zeros(Ny+1,Nx);

% Wash Channel
for i = 2:Ny+1
for j = 2:NxL
    dy = YV(i,j) - YV(i-1,j);

    % Top and Bottom Edges of cell
    ys = YV(i-1,j); yn = YV(i,j);

    % Location of Nodes
    XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
    YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);


    % East
    ue = 0.5*(U(i,j+1) + U(i,j));
    me = 0.5 * dy * (U(i,j) + U(i,j+1));

    % North
    if i == Ny+1 
        un = U(i+1,j); % Don't interpolate at boundary
    else
        un = U(i,j) + (yn - YP) * ( U(i+1,j) - U(i,j) ) / ( YN - YP );
    end
    mn = 0.5 * ( V(i,j) *(XP - XW) + V(i,j+1) *(XE - XP) );


    % West
    uw = 0.5*(U(i,j) + U(i,j-1));
    mw = 0.5 * dy * (U(i,j) + U(i,j-1));

    % South
    if i == 2
        us = U(i-1,j); % Don't interpolate at boundary
    else
        us = U(i,j) + (ys-YP) * ( U(i-1,j) - U(i,j) ) / (YS-YP);
    end
    ms = 0.5 * ( V(i-1,j)*(XP-XW) + V(i-1,j+1)*(XE-XP) );

    NLU(i,j) = me*ue +mn*un -mw*uw - ms*us;
end
end

% Injury Channel
for i = NyL+2:NyM+1
for j = NxL+1:NxM+1
    dy = YV(i,j) - YV(i-1,j);

    % Top and Bottom Edges of cell
    ys = YV(i-1,j); yn = YV(i,j);

    % Location of Nodes
    XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
    YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);


    % East
    ue = 0.5*(U(i,j+1) + U(i,j));
    me = 0.5 * dy * (U(i,j) + U(i,j+1));

    % North
    if (i == NyM+1) && (j ~= NxL+1)
        un = U(i+1,j); % Don't interpolate at boundary
    elseif (i == NyM+1) && (j ~= NxM+1)
        un = U(i+1,j); % Don't interpolate at boundary
    else
        un = U(i,j) + (yn - YP) * ( U(i+1,j) - U(i,j) ) / ( YN - YP );
    end
    mn = 0.5 * ( V(i,j) *(XP - XW) + V(i,j+1) *(XE - XP) );


    % West
    uw = 0.5*(U(i,j) + U(i,j-1));
    mw = 0.5 * dy * (U(i,j) + U(i,j-1));

    % South
    if i == NyL+2
        us = U(i-1,j); % Don't interpolate at boundary
    else
        us = U(i,j) + (ys-YP) * ( U(i-1,j) - U(i,j) ) / (YS-YP);
    end
    ms = 0.5 * ( V(i-1,j)*(XP-XW) + V(i-1,j+1)*(XE-XP) );

    NLU(i,j) = me*ue +mn*un -mw*uw - ms*us;
end
end

% Blood Channel
for i = 2:Ny+1
for j = NxM+2:Nx
    dy = YV(i,j) - YV(i-1,j);

    % Top and Bottom Edges of cell
    ys = YV(i-1,j); yn = YV(i,j);

    % Location of Nodes
    XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
    YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);


    % East
    ue = 0.5*(U(i,j+1) + U(i,j));
    me = 0.5 * dy * (U(i,j) + U(i,j+1));

    % North
    if i == Ny+1 
        un = U(i+1,j); % Don't interpolate at boundary
    else
        un = U(i,j) + (yn - YP) * ( U(i+1,j) - U(i,j) ) / ( YN - YP );
    end
    mn = 0.5 * ( V(i,j) *(XP - XW) + V(i,j+1) *(XE - XP) );


    % West
    uw = 0.5*(U(i,j) + U(i,j-1));
    mw = 0.5 * dy * (U(i,j) + U(i,j-1));

    % South
    if i == 2
        us = U(i-1,j); % Don't interpolate at boundary
    else
        us = U(i,j) + (ys-YP) * ( U(i-1,j) - U(i,j) ) / (YS-YP);
    end
    ms = 0.5 * ( V(i-1,j)*(XP-XW) + V(i-1,j+1)*(XE-XP) );

    NLU(i,j) = me*ue +mn*un -mw*uw - ms*us;
end
end

% --- Compute NLV for y-momentum equation --------------------------------
NLV = zeros(Ny,Nx+1);

% Wash Channel
for i = 2:Ny
for j = 2:NxL+1

    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);

    % Left and Right Edges of cell
    xw = XU(i,j-1); xe = XU(i,j);

    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);

    % East
    if (j == NxL+1) && (i < NyL+1)
        ve = V(i,j+1); % don't interpolate on boundary
    elseif (j == NxL+1) && (i > NyM+1)
        ve = V(i,j+1); % don't interpolate on boundary
    else
        ve = V(i,j) + (xe-XP)* (V(i,j+1) - V(i,j)) / (XE-XP);
    end
    me = 0.5*( U(i+1,j) * (YN-YP) + U(i,j) * (YP-YS) );

    % North
    vn = 0.5*( V(i+1,j) + V(i,j) );
    mn = 0.5* dx * ( V(i+1,j) + V(i,j) );


    % West
    if j == 2
        vw = V(i,j-1); % don't interpolate on boundary
    else
        vw = V(i,j) + (xw-XP) * ( V(i,j-1)-V(i,j) ) / (XW-XP);
    end
    mw = 0.5*( U(i+1,j-1) * (YN-YP) + U(i,j-1) * (YP-YS) );

    % South
    vs = 0.5*( V(i,j) + V(i-1,j) );
    ms = 0.5*dx * ( V(i,j) + V(i-1,j) );

    NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;
end
end

% Injury Channel
for i = NyL+2:NyM
for j = NxL+2:NxM+1

    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);

    % Left and Right Edges of cell
    xw = XU(i,j-1); xe = XU(i,j);

    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);

    % East
    ve = V(i,j) + (xe-XP)* (V(i,j+1) - V(i,j)) / (XE-XP);
    me = 0.5*( U(i+1,j) * (YN-YP) + U(i,j) * (YP-YS) );

    % North
    vn = 0.5*( V(i+1,j) + V(i,j) );
    mn = 0.5* dx * ( V(i+1,j) + V(i,j) );


    % West
    vw = V(i,j) + (xw-XP) * ( V(i,j-1)-V(i,j) ) / (XW-XP);
    mw = 0.5*( U(i+1,j-1) * (YN-YP) + U(i,j-1) * (YP-YS) );

    % South
    vs = 0.5*( V(i,j) + V(i-1,j) );
    ms = 0.5*dx * ( V(i,j) + V(i-1,j) );

    NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;
end
end

% Wash Channel
for i = 2:Ny
for j = NxM+2:Nx+1

    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);

    % Left and Right Edges of cell
    xw = XU(i,j-1); xe = XU(i,j);

    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);

    % East
    if j == Nx+1
        ve = V(i,j+1); % don't interpolate on boundary
    else
        ve = V(i,j) + (xe-XP)* (V(i,j+1) - V(i,j)) / (XE-XP);
    end
    me = 0.5*( U(i+1,j) * (YN-YP) + U(i,j) * (YP-YS) );

    % North
    vn = 0.5*( V(i+1,j) + V(i,j) );
    mn = 0.5* dx * ( V(i+1,j) + V(i,j) );


    % West
    if (j == NxM+2) && (i < NyL+1) 
        vw = V(i,j-1); % don't interpolate on boundary
    elseif (j == NxM+2) && (i > NyM+1)
        vw = V(i,j-1); % don't interpolate on boundary
    else
        vw = V(i,j) + (xw-XP) * ( V(i,j-1)-V(i,j) ) / (XW-XP);
    end
    mw = 0.5*( U(i+1,j-1) * (YN-YP) + U(i,j-1) * (YP-YS) );

    % South
    vs = 0.5*( V(i,j) + V(i-1,j) );
    ms = 0.5*dx * ( V(i,j) + V(i-1,j) );

    NLV(i,j) = me*ve + mn*vn - mw*vw - ms*vs;
end
end
