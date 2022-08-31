function [A] = BuildAuH(N,XU,YU,XV,YV,Re,dt)
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
% Au is NNXY x NNXY matrix
NNXY = (Nx+1)*(Ny+2) - (NxM - NxL -1) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 1) + (Nx + 1 - NxM); 

% Initialize the counter used in building the sparse matrices
ctr = 0;

% row = zeros(5*NNXY,1); col = row; A = row;

% --- Boundary Conditions ------------------------------------------------

% Top Boundaries (\partial \Omega_i, i = 1, 3)
for j = 1:Nsml
    ctr = ctr + 1;
    row(ctr) = NNXY - Nsml + j;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end

% Bottom Left Boundary (\partial \Omega_2)
for j = 1:NxL+1
    % Principal - Boundary
    ctr = ctr + 1;
    row(ctr) = j;
    col(ctr) = row(ctr);
    A(ctr) = -1/(YU(2,j)-YU(1,j));
    
    % North
    ctr = ctr + 1;
    row(ctr) = j;
    col(ctr) = row(ctr) + Nsml;
    A(ctr) = 1/(YU(2,j)-YU(1,j));
end

% Bottom Right Boundary (\partial \Omega_4)
j_ctr = NxL+1;
for j = NxM+1:Nx+1
    j_ctr = j_ctr+1;
    
    % Principal - Boundary
    ctr = ctr + 1;
    row(ctr) = j_ctr;
    col(ctr) = row(ctr);
    A(ctr) = -1/(YU(2,j)-YU(1,j));
    
    % North
    ctr = ctr + 1;
    row(ctr) = j_ctr;
    col(ctr) = row(ctr) + Nsml;
    A(ctr) = 1/(YU(2,j)-YU(1,j));
end

% Left and Right Boundaries (\partial \Omega_i, i = 5,6)
% a_i * Up + b_i * (dUp/dt + Uinf * dUp/dx) = g_i
for i = 2:Ny+1
    if i < NyL + 1 % Lower Legs
        ROW_L = (i-1)*Nsml + 1;
        ROW_R = (i)*Nsml;
    elseif (i > NyL) && (i < NyM + 3) % Middle of H
        ROW_L = NyL*Nsml + (i - NyL - 1) * (Nx+1) + 1;
        ROW_R = NyL*Nsml + (i - NyL) * (Nx+1);
    else % Upper Legs
        ROW_L = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (i-NyM - 3)*Nsml + 1;
        ROW_R = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (i-NyM - 2)*Nsml;
    end
    
    % Left \partial \Omega_5
    ctr = ctr + 1;
    row(ctr) = ROW_L;
    col(ctr) = row(ctr);
    A(ctr) = 1;
    
    % Right \partial \Omega_6
    ctr = ctr + 1;
    row(ctr) = ROW_R;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end

% Left and Right Walls for Upper Legs (\Gamma_i, i = 1,3)
for i = NyM + 2 : Ny + 1
    if i == NyM + 2
        ROW_L = NyL*Nsml + (NyM+1 - NyL)*(Nx+1) + NxL + 1;
        ROW_R = NyL*Nsml + (NyM+1 - NyL)*(Nx+1) + NxM + 1;
    else
        ROW_L = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) ...
                + (i-NyM - 3)*Nsml + NxL + 1;
        ROW_R = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) ...
                + (i-NyM - 3)*Nsml + NxL + 2;
    end

    % Left \Gamma_1
    ctr = ctr + 1;
    row(ctr) = ROW_L;
    col(ctr) = row(ctr);
    A(ctr) = 1;
    
    % Right \Gamma_3
    ctr = ctr + 1;
    row(ctr) = ROW_R;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end

% Left and Right Walls for Lower Legs (\Gamma_i, i = 2,4)
for i = 2:NyL+1
    if i < NyL + 1
        ROW_L = (i-1)*Nsml + NxL + 1;
        ROW_R = (i-1)*Nsml + NxL + 2;
    else
        ROW_L = (i-1)*Nsml + NxL + 1;
        ROW_R = (i-1)*Nsml + NxM + 1;
    end

    % Left \Gamma_2
    ctr = ctr + 1;
    row(ctr) = ROW_L;
    col(ctr) = row(ctr);
    A(ctr) = 1;
    
    % Right \Gamma_4
    ctr = ctr + 1;
    row(ctr) = ROW_R;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end


% Upper Mid-H Wall \Gamma_5
for j = NxL + 2 : NxM
    % North - Ghost
    ctr = ctr + 1;
    row(ctr) = NyL*Nsml + (NyM + 1 - NyL) * (Nx+1) + j;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end

% Lower Mid-H Wall \Gamma_6
for j = NxL + 2 : NxM
    % South - Ghost
    ctr = ctr + 1;
    row(ctr) = NyL*Nsml + j;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end

% --- End of Boundary Conditions ----------------------------------------



% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    j_ctr = NxL+2; % Starting value for index on right of H
    for j = [2:NxL, NxM+2:Nx]
        dy = YV(i,j)-YV(i-1,j);
        yn = YU(i+1,j);
        yp = YU(i,j);
        ys = YU(i-1,j);
        dx = XV(i,j+1) - XV(i,j);
        xe = XU(i,j+1);
        xp = XU(i,j);
        xw = XU(i,j-1);
        
        if j < NxM+2
            ROW = Nsml*(i-1)+j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        end
        
        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) + Nsml;
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) - Nsml;
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr);
        A(ctr) = -DE -DN -DW -DS + dx*dy/dt;
    
    end
end

% ---- Fill Inner cells at Lower Interface yu_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
j_ctr = NxL+2; % Starting matrix value right H when i == NyL only 

for j = [2:NxL, NxM+2:Nx]
    dy = YV(i,j)-YV(i-1,j);
    yn = YU(i+1,j);
    yp = YU(i,j);
    ys = YU(i-1,j);
    dx = XV(i,j+1) - XV(i,j);
    xe = XU(i,j+1);
    xp = XU(i,j);
    xw = XU(i,j-1);
    if (i == NyL) && j > NxM+1
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    else
        ROW = Nsml*(i-1)+j;
    end

    % East
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr) + 1;
    DE = - 1/Re * dy/(xe-xp)/2; 
    A(ctr) = DE;

    % North
    ctr = ctr + 1;
    row(ctr) = ROW;
    if (i == NyL) && (j < NxL + 1)
        col(ctr) = row(ctr) + Nsml;
    else 
        col(ctr) = row(ctr) + (Nx+1);
    end
    DN = - 1/Re * dx/(yn-yp)/2;  
    A(ctr) = DN;

    % West
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr) - 1;
    DW = - 1/Re * dy/(xp-xw)/2;
    A(ctr) = DW;

    % South
    ctr = ctr + 1;
    row(ctr) = ROW;
    if i == NyL
        col(ctr) = row(ctr) - Nsml;
    elseif (i == NyL+1) && (j < NxL + 1)
        col(ctr) = row(ctr) - Nsml;
    else
        col(ctr) = row(ctr) - (Nx+1);
    end
    DS = - 1/Re * dx/(yp-ys)/2;
    A(ctr) = DS;

    % Principal
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr);
    A(ctr) = -DE -DN -DW -DS + dx*dy/dt;

end
end

% ---- Fill inner cells in middle of H for i = NyL +2, ..., NyM+1
for i = NyL+2 : NyM+1
for j = 2:Nx
    dy = YV(i,j)-YV(i-1,j);
    yn = YU(i+1,j);
    yp = YU(i,j);
    ys = YU(i-1,j);
    dx = XV(i,j+1) - XV(i,j);
    xe = XU(i,j+1);
    xp = XU(i,j);
    xw = XU(i,j-1);
    
    ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;

    % East
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr) + 1;
    DE = - 1/Re * dy/(xe-xp)/2; 
    A(ctr) = DE;

    % North
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr) + (Nx+1);
    DN = - 1/Re * dx/(yn-yp)/2;  
    A(ctr) = DN;

    % West
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr) - 1;
    DW = - 1/Re * dy/(xp-xw)/2;
    A(ctr) = DW;

    % South
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr) - (Nx+1);
    DS = - 1/Re * dx/(yp-ys)/2;
    A(ctr) = DS;

    % Principal
    ctr = ctr + 1;
    row(ctr) = ROW;
    col(ctr) = row(ctr);
    A(ctr) = -DE -DN -DW -DS + dx*dy/dt;
end
end

% ---- Fill Inner cells at Upper Interface yu_i, i = NyM+2, NyM+3 ----
for i = NyM+2 : NyM+3

j_ctr = NxL+2; % Starting matrix value right H when i == NyM+3 only 
    for j = [2:NxL, NxM+2:Nx]
        dy = YV(i,j)-YV(i-1,j);
        yn = YU(i+1,j);
        yp = YU(i,j);
        ys = YU(i-1,j);
        dx = XV(i,j+1) - XV(i,j);
        xe = XU(i,j+1);
        xp = XU(i,j);
        xw = XU(i,j-1);
        
        if (i == NyM+3) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + j_ctr;
        else
            ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;
        end

        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;

        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        if (i == NyM+2) && (j < NxL + 1)
            col(ctr) = row(ctr) + (Nx+1);    
        else 
            col(ctr) = row(ctr) + Nsml;
        end
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;

        % South
        ctr = ctr + 1;
        row(ctr) = ROW;
        if i == NyM+3 && (j > NxM+1)
            col(ctr) = row(ctr) - Nsml;
        else
            col(ctr) = row(ctr) - (Nx+1);
        end
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr);
        A(ctr) = -DE -DN -DW -DS + dx*dy/dt;

    end
end

% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+4 : Ny+1
    j_ctr = NxL+2; % Starting value for index on right of H
    for j = [2:NxL, NxM+2:Nx]
        dy = YV(i,j)-YV(i-1,j);
        yn = YU(i+1,j);
        yp = YU(i,j);
        ys = YU(i-1,j);
        dx = XV(i,j+1) - XV(i,j);
        xe = XU(i,j+1);
        xp = XU(i,j);
        xw = XU(i,j-1);
        
        if j < NxM+2
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j_ctr;
        end
        
        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) + Nsml;
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr) - Nsml;
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = row(ctr);
        A(ctr) = -DE -DN -DW -DS +  dx*dy/dt;
    end
end
% -----------------------------------------------------------------------

% Remove zero entries from row, col, A
% A = nonzeros(A); 
% row = nonzeros(row); 
% col = nonzeros(col);

% Build Sparse Matrix
A = sparse(row, col, A, NNXY, NNXY);