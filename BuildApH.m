function[Ap] = BuildApH(N,xp,yp,XU,YV)

% Note: i,j -> Grid locations while i_ctr,j_ctr -> Matrix Locations

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Ap is NNXY x NNXY matrix without Augmentation
NNXY = (Nx)*(Ny) - (NxM - NxL) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsml x Nsml
Nsml = (NxL) + (Nx - NxM); 

% Initialize the counter used in building the sparse matrices
ctr = 0;

% ---- Fill Cells in Lower Legs of H ------------------------------
for i = 1:NyL
    
    j_ctr = NxL; % Starting value for index when j > NxL
    
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*(i-1)+j_ctr;
    end
    
    % Location of grid points (Note xp & yp have locations of ghost nodes)
    xj = xp(j+1);
    xe = xp(j+2);
    xw = xp(j);
    yi = yp(i+1);
    yn = yp(i+2);
    ys = yp(i);
    dx = XU(i+1,j+1)-XU(i+1,j);
    dy = YV(i+1,j+1)-YV(i,j+1);

    % Coefficients
    ae = 1/(xe-xj)/dx;
    aw = 1/(xj-xw)/dx;
    an = 1/(yn-yi)/dy;
    as = 1/(yi-ys)/dy; 
    Bp = 0; Lp = 0; Rp = 0;
    
    % Extra Coefficients for Ghost Nodes
    if i == 1 % On \partial \Omega_2 or \partial \Omega_4
        Bp = -as;
    end 
    if (j == 1) || (j == NxM+1) % On \partial \Omega_5 or \Gamma_4 boundary
        Lp = aw;
    end 
    if (j == NxL) || (j == Nx) % On \Gamma_2 or \partial \Omega_6 boundary
        Rp = ae;
    end
    
    % Principal
    ctr = ctr+1;
    ap =  -(ae + aw  + an + as) + Bp + Lp + Rp;
    row(ctr) = ROW_ind;
    col(ctr) = row(ctr);
    A(ctr)   = ap;
    
    % East - Wash Channel
    if j < NxL
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = ae;
    end
    
    % East - Blood Channel
    if (j > NxM) && (j < Nx)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = ae;
    end

    % West - Wash Channel
    if (j > 1) && (j < NxL+1)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = aw;
    end
    
    % West - Blood Channel
    if (j > NxM+1) 
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = aw;
    end
    
    % North
    if (i == NyL) && (j > NxL) % At interface of lower legs & injury chan.
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nx;
        A(ctr)   = an;
    else
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nsml;
        A(ctr)   = an;
    end
    
    % South
    if i > 1
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nsml;
        A(ctr)   = as;
    end
   
end
end

% ---- Fill Cells in Injury Channel ---
for i = NyL+1 : NyM   
for j = 1:Nx
    ROW_ind = Nsml*(NyL) + (i - NyL -1) * Nx + j;
    
    % Location of grid points (Note xp & yp have locations of ghost nodes)
    xj = xp(j+1);
    xe = xp(j+2);
    xw = xp(j);
    yi = yp(i+1);
    yn = yp(i+2);
    ys = yp(i);
    dx = XU(i+1,j+1)-XU(i+1,j);
    dy = YV(i+1,j+1)-YV(i,j+1);

    % Coefficients
    ae = 1/(xe-xj)/dx;
    aw = 1/(xj-xw)/dx;
    an = 1/(yn-yi)/dy;
    as = 1/(yi-ys)/dy; 
    Bp = 0; Tp = 0; Lp = 0; Rp = 0;
    
    % Extra Coefficients for Ghost Nodes
    if i == (NyL+1) 
        if (j > NxL) && (j < NxM+1) % On \Gamma_6 boundary
            Bp = as;
        end
    end 
    if (j == 1) % On \partial \Omega_5 
        Lp = aw;
    end 
    if (j == Nx) % On \partial \Omega_6 boundary
        Rp = ae;
    end
    if i == (NyM) 
        if (j > NxL) && (j < NxM+1) % On \Gamma_5 boundary
            Tp = an;
        end
    end 
    
    % Principal
    ctr = ctr+1;
    ap =  -(ae + aw  + an + as) + Bp + Tp + Lp + Rp;
    row(ctr) = ROW_ind;
    col(ctr) = row(ctr);
    A(ctr)   = ap;
    
    % East 
    if j < Nx
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = ae;
    end
    
    % West 
    if (j > 1) 
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = aw;
    end
    
    % North
    if (i < NyM) % in injury channel
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nx;
        A(ctr)   = an;
    elseif (i == NyM) && (j < NxL+1) % at interface left
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nx;
        A(ctr)   = an;
    elseif (i == NyM) && (j > NxM) % at interface right
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nsml;
        A(ctr)   = an;
    end
    
    % South
    if (i > NyL+1) % in injury channel
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nx;
        A(ctr)   = as;
    elseif (i == NyL+1) && (j < NxL+1)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nsml;
        A(ctr)   = as;
    elseif (i == NyL+1) && (j > NxM)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nx;
        A(ctr)   = as;
    end
 
end
end

% ---- Fill Cells in Upper Legs of H ------------------------------
for i = NyM+1:Ny
    
    j_ctr = NxL; % Starting value for index when j > NxL
    
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j_ctr;
    end
    
    % Location of grid points (Note xp & yp have locations of ghost nodes)
    xj = xp(j+1);
    xe = xp(j+2);
    xw = xp(j);
    yi = yp(i+1);
    yn = yp(i+2);
    ys = yp(i);
    dx = XU(i+1,j+1)-XU(i+1,j);
    dy = YV(i+1,j+1)-YV(i,j+1);

    % Coefficients
    ae = 1/(xe-xj)/dx;
    aw = 1/(xj-xw)/dx;
    an = 1/(yn-yi)/dy;
    as = 1/(yi-ys)/dy; 
    Tp = 0; Lp = 0; Rp = 0;
    
    % Extra Coefficients for Ghost Nodes
    if (i == Ny) % On \partial \Omega_1 or \partial \Omega_3
        Tp = -an;
    end 
    if (j == 1) || (j == NxM+1) % On \partial \Omega_5 or \Gamma_3 boundary
        Lp = aw;
    end 
    if (j == NxL) || (j == Nx) % On \Gamma_1 or \partial \Omega_6 boundary
        Rp = ae;
    end
    
    % Principal
    ctr = ctr+1;
    ap =  -(ae + aw  + an + as) + Tp + Lp + Rp;
    row(ctr) = ROW_ind;
    col(ctr) = row(ctr);
    A(ctr)   = ap;
    
    % East - Wash Channel
    if j < NxL
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = ae;
    end
    
    % East - Blood Channel
    if (j > NxM) && (j < Nx)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + 1;
        A(ctr)   = ae;
    end

    % West - Wash Channel
    if (j > 1) && (j < NxL+1)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = aw;
    end
    
    % West - Blood Channel
    if (j > NxM+1) 
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - 1;
        A(ctr)   = aw;
    end
    
    % North
    if (i < Ny)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) + Nsml;
        A(ctr)   = an;
    end
    
    % South
    if (i == NyM+1) && (j < NxL+1)
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nx;
        A(ctr)   = as;
    else
        ctr = ctr+1;
        row(ctr) = ROW_ind;
        col(ctr) = row(ctr) - Nsml;
        A(ctr)   = as;
    end
   
end
end

Ap = sparse(row,col,A);
