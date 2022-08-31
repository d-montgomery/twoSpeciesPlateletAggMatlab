function [WU,WV] = hinderedUV_UpW(N,U,V,W)

% Calculate the hindered velocity field by interpolating W onto U and V

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Initialize hindered velocity fields WU and WV
WU = U; 
WV = V;

% --- Compute WU ---

% Wash channel
for i = 2:Ny+1
for j = 2:NxL
    WU(i,j) = max(0,U(i,j))*W(i,j) + min(0,U(i,j))*W(i,j+1);
end
end

% Injury Channel
for i = NyL+2:NyM+1
for j = NxL+1:NxM+1
    WU(i,j) = max(0,U(i,j))*W(i,j) + min(0,U(i,j))*W(i,j+1);
end
end

% Blood channel
for i = 2:Ny+1
for j = NxM+2:Nx
    WU(i,j) = max(0,U(i,j))*W(i,j) + min(0,U(i,j))*W(i,j+1);
end
end

% --- Compute WV ---

% Wash channel
for i = 2:Ny
for j = 2:NxL+1
    WV(i,j) = max(0,V(i,j))*W(i,j) + min(0,V(i,j))*W(i+1,j);
end
end

% Injury channel
for i = NyL+2:NyM
for j = NxL+2:NxM+1
    WV(i,j) = max(0,V(i,j))*W(i,j) + min(0,V(i,j))*W(i+1,j);
end
end

% Blood channel
for i = 2:Ny
for j = NxM+2:Nx+1
    WV(i,j) = max(0,V(i,j))*W(i,j) + min(0,V(i,j))*W(i+1,j);
end
end
