function [A] = BuildBvH(N,XU,YU,YV,B)
% INPUTS:
% Nx and Ny are the total number of cells in x-y directions
% Nxs and Nys are the number of cells in the step
% xu and yu is where v(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% a and b are the BC coefficients
% Re is the Reynolds number
% dt is time step

% Note: i,j -> Grid locations while i_ctr,j_ctr -> Matrix Locations

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Av is NNXY x NNXY matrix
NNXY = (Nx+2)*(Ny+1) - (NxM - NxL -2) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 2) + (Nx + 2 - NxM); 

% Initialize the counter used in building the sparse matrices
ctr = 0;


% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    
    j_ctr = NxL+3; % Starting value for index on right of H
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        [Bp,dx,dy] = Bp_Vcoeff(i,j,XU,YU,YV,B);
        
        if j < NxM+2
            ROW = Nsml*(i-1)+j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        end

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = dx*dy*Bp/2;
    
    end
end


% ---- Fill Inner cells at Lower Interface yv_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
  
    j_ctr = NxL+3; % Starting matrix value right H when i == NyL only
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        [Bp,dx,dy] = Bp_Vcoeff(i,j,XU,YU,YV,B);
        
        if (i == NyL) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        else
            ROW = Nsml*(i-1)+j;
        end
        
        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = dx*dy*Bp/2;
    
    end
end


% ---- Fill inner cells in middle of H for i = NyL +2, ..., NyM
for i = NyL+2 : NyM
    for j = 2 : Nx+1
        [Bp,dx,dy] = Bp_Vcoeff(i,j,XU,YU,YV,B);
        
        ROW = NyL*Nsml + (i - NyL-1)*(Nx+2) + j;
        
        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = dx*dy*Bp/2;
    
    end
end


% ---- Fill Inner cells at Upper Interface yv_i, i = NyM+1, NyM+2 ---
for i = NyM+1 : NyM+2
    j_ctr = NxL+3; % Starting matrix value right H when i == NyM+2 only
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        [Bp,dx,dy] = Bp_Vcoeff(i,j,XU,YU,YV,B);
        
        if (i == NyM+2) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+1 - NyL) * (Nx+2) + j_ctr;
        else
            ROW = Nsml*NyL + (i-1 - NyL) * (Nx+2)+j;
        end
        
        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = dx*dy*Bp/2;
    
    end
end

% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+3 : Ny

    j_ctr = NxL+3; % Starting value for index on right of H
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        [Bp,dx,dy] = Bp_Vcoeff(i,j,XU,YU,YV,B);
        
        if j < NxM+2
            ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j_ctr;
        end
        
        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = dx*dy*Bp/2;
    
    end
end

A = sparse(row, col, A, NNXY, NNXY);