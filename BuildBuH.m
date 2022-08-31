function [A] = BuildBuH(N,XU,XV,YV,B)
% INPUTS:
% N = [NxL, NxM, Nx, NyL, NyM, Ny];
% xu and yu is where v(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid


% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Au is NNXY x NNXY matrix
NNXY = (Nx+1)*(Ny+2) - (NxM - NxL -1) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 1) + (Nx + 1 - NxM); 

% Initialize the counter used in building the sparse matrices
ctr = 0;

% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    j_ctr = NxL+2; % Starting value for index on right of H
    for j = [2:NxL, NxM+2:Nx]
        [Bp,dx,dy] = Bp_Ucoeff(i,j,XU,XV,YV,B);
        
        if j < NxM+2
            ROW = Nsml*(i-1)+j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        end
        

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr);
        A(ctr) = dx*dy*Bp/2;
    
    end
end

% ---- Fill Inner cells at Lower Interface yu_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
j_ctr = NxL+2; % Starting matrix value right H when i == NyL only 

for j = [2:NxL, NxM+2:Nx]
    [Bp,dx,dy] = Bp_Ucoeff(i,j,XU,XV,YV,B);
    
    if (i == NyL) && j > NxM+1
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    else
        ROW = Nsml*(i-1)+j;
    end

    % Principal
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr);
    A(ctr) = dx*dy*Bp/2;

end
end

% ---- Fill inner cells in middle of H for i = NyL +2, ..., NyM+1
for i = NyL+2 : NyM+1
for j = 2:Nx
    [Bp,dx,dy] = Bp_Ucoeff(i,j,XU,XV,YV,B);
    
    ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;

    % Principal
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr);
    A(ctr) = dx*dy*Bp/2;
end
end

% ---- Fill Inner cells at Upper Interface yu_i, i = NyM+2, NyM+3 ----
for i = NyM+2 : NyM+3

j_ctr = NxL+2; % Starting matrix value right H when i == NyM+3 only 
    for j = [2:NxL, NxM+2:Nx]
        [Bp,dx,dy] = Bp_Ucoeff(i,j,XU,XV,YV,B);
        
        if (i == NyM+3) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + j_ctr;
        else
            ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;
        end

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr);
        A(ctr) = dx*dy*Bp/2;

    end
end

% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+4 : Ny+1
    j_ctr = NxL+2; % Starting value for index on right of H
    for j = [2:NxL, NxM+2:Nx]
        [Bp,dx,dy] = Bp_Ucoeff(i,j,XU,XV,YV,B);
        
        if j < NxM+2
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j_ctr;
        end

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr);
        A(ctr) = dx*dy*Bp/2;
    end
end
% -----------------------------------------------------------------------

% Build Sparse Matrix
A = sparse(row, col, A, NNXY, NNXY);

end 
