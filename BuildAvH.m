function [A] = BuildAvH(N,XU,YU,XV,YV,Re,dt)
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

row = zeros(5*NNXY,1); col = row; A = row;

% --- Boundary Conditions ------------------------------------------------

% Top Boundary (\partial \Omega_i, i = 1, 3)
% c_i * Vp + d_i * (dVp/dt + Vinf * dVp/dy) = h_i
j_ctr = NxL+2; % indexing for right bc starts at NxL + 4
for j = [1 : NxL+2, NxM+1 : Nx+2]
    if j < NxL + 3 % top left bc
        % Principal
        ctr = ctr + 1;
        row(ctr) = NNXY - Nsml + j;
        col(ctr) = row(ctr);
        A(ctr) = 1;
    else % top right bc
        % Principal
        j_ctr = j_ctr +1;
        ctr = ctr + 1;
        row(ctr) = NNXY - Nsml + j_ctr;
        col(ctr) = row(ctr);
        A(ctr) = 1;
    end
end

% Bottom Boundaries (\partial \Omega_i, i = 2, 4)
% c_i * Vp + d_i * (dVp/dt + Vinf * dVp/dy) = h_i
j_ctr = NxL+2; % indexing for right bc starts at NxL + 4
for j = [1 : NxL+2, NxM+1 : Nx+2]
    if j < NxL + 3 % bottom left bc
        % Principal
        ctr = ctr + 1;
        row(ctr) = j;
        col(ctr) = row(ctr);
        A(ctr) = -1/( YV(2,j) - YV(1,j) );
        
        % North
        ctr = ctr + 1;
        row(ctr) = j;
        col(ctr) = row(ctr) + Nsml;
        A(ctr) = 1/( YV(2,j) - YV(1,j) );
        
    else % bottom right bc
        % Principal
        j_ctr = j_ctr +1;
        ctr = ctr + 1;
        row(ctr) = j_ctr;
        col(ctr) = row(ctr);
        A(ctr) = - 1/( YV(2,j) - YV(1,j) );
        
        % North
        ctr = ctr + 1;
        row(ctr) = j_ctr;
        col(ctr) = row(ctr)+Nsml;
        A(ctr) = 1/( YV(2,j) - YV(1,j) );

    end
end



% Left and Right Boundaries (\partial \Omega_i, i = 5,6)
for i = 2:Ny
    if i < NyL + 1 % Lower Legs i = 1:NyL
        ROW_L = (i-1)*Nsml + 1;
        ROW_R = (i)*Nsml;
    elseif (i > NyL) && (i < NyM + 2) % Middle of H i = NyL+1:NyM+1
        ROW_L = NyL*Nsml + (i - NyL - 1) * (Nx+2) + 1;
        ROW_R = NyL*Nsml + (i - NyL) * (Nx+2);
    else % Upper Legs i = NyM+2:Ny+1
        ROW_L = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (i-NyM - 2)*Nsml + 1;
        ROW_R = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (i-NyM - 1)*Nsml;
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
for i = NyM + 2 : Ny 
        ROW_L = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) ...
                + (i-NyM - 2)*Nsml + NxL + 2;
        ROW_R = ROW_L + 1;

    % \Gamma_1 - Upper Left
    % East - Ghost
    ctr = ctr + 1;
    row(ctr) = ROW_L;
    col(ctr) = row(ctr);
    A(ctr) = 1;
    
    
    % \Gamma_3 - Upper Right
    % West - Ghost
    ctr = ctr + 1;
    row(ctr) = ROW_R;
    col(ctr) = row(ctr);
    A(ctr) = 1;

end

% Left and Right Walls for Lower Legs (\Gamma_i, i = 2,4)
for i = 2:NyL 
    ROW_L = (i-1)*Nsml + NxL + 2;
    ROW_R = ROW_L + 1;

    % \Gamma_2 - Lower Left
    % East - Ghost
    ctr = ctr + 1;
    row(ctr) = ROW_L;
    col(ctr) = row(ctr);
    A(ctr) = 1;
 
    
    % \Gamma_4 - Lower Right
    % West - Ghost
    ctr = ctr + 1;
    row(ctr) = ROW_R;
    col(ctr) = row(ctr);
    A(ctr) = 1;

end

% Upper and lower Mid-H Wall \Gamma_i, i = 5, 6
for j = NxL + 2 : NxM + 1
    
    % Upper wall \Gamma 5
    ctr = ctr + 1;
    row(ctr) = NyL*Nsml + (NyM - NyL) * (Nx+2) + j;
    col(ctr) = row(ctr);
    A(ctr) = 1;
    
    % Upper wall \Gamma 6
    ctr = ctr + 1;
    row(ctr) = NyL*Nsml + j;
    col(ctr) = row(ctr);
    A(ctr) = 1;
end

% --- End of Boundary Conditions ----------------------------------------



% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    
    j_ctr = NxL+3; % Starting value for index on right of H
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        dy = YU(i+1,j)-YU(i,j);
        yn = YV(i+1,j);
        yp = YV(i,j);
        ys = YV(i-1,j);
        dx = XU(i,j) - XU(i,j-1);
        xe = XV(i,j+1);
        xp = XV(i,j);
        xw = XV(i,j-1);
        
        if j < NxM+2
            ROW = Nsml*(i-1)+j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        end
        
        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + Nsml;
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - Nsml;
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = -DE -DN -DW -DS + dx*dy/dt;
    
    end
end


% ---- Fill Inner cells at Lower Interface yv_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
  
    j_ctr = NxL+3; % Starting matrix value right H when i == NyL only
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        dy = YU(i+1,j)-YU(i,j);
        yn = YV(i+1,j);
        yp = YV(i,j);
        ys = YV(i-1,j);
        dx = XU(i,j) - XU(i,j-1);
        xe = XV(i,j+1);
        xp = XV(i,j);
        xw = XV(i,j-1);
        
        if (i == NyL) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        else
            ROW = Nsml*(i-1)+j;
        end
        
        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        if (i == NyL) && (j < NxL + 2)
            col(ctr) = row(ctr) + Nsml;
        else 
            col(ctr) = row(ctr) + (Nx+2);
        end
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South 
        ctr = ctr + 1;
        row(ctr) = ROW;
        if i == NyL
            col(ctr) = row(ctr) - Nsml;
        elseif (i == NyL+1) && (j < NxL + 2)
            col(ctr) = row(ctr) - Nsml;
        else
                col(ctr) = row(ctr) - (Nx+2);
        end
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = -DE -DN -DW -DS + dx*dy/dt;
    
    end
end


% ---- Fill inner cells in middle of H for i = NyL +2, ..., NyM
for i = NyL+2 : NyM
    for j = 2 : Nx+1
        dy = YU(i+1,j)-YU(i,j);
        yn = YV(i+1,j);
        yp = YV(i,j);
        ys = YV(i-1,j);
        dx = XU(i,j) - XU(i,j-1);
        xe = XV(i,j+1);
        xp = XV(i,j);
        xw = XV(i,j-1);
        
        ROW = NyL*Nsml + (i - NyL-1)*(Nx+2) + j;
        
        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + (Nx+2);
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - (Nx+2);
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = -DE -DN -DW -DS + dx*dy/dt;
    
    end
end


% ---- Fill Inner cells at Upper Interface yv_i, i = NyM+1, NyM+2 ---
for i = NyM+1 : NyM+2
    j_ctr = NxL+3; % Starting matrix value right H when i == NyM+2 only
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        dy = YU(i+1,j)-YU(i,j);
        yn = YV(i+1,j);
        yp = YV(i,j);
        ys = YV(i-1,j);
        dx = XU(i,j) - XU(i,j-1);
        xe = XV(i,j+1);
        xp = XV(i,j);
        xw = XV(i,j-1);
        
        if (i == NyM+2) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+1 - NyL) * (Nx+2) + j_ctr;
        else
            ROW = Nsml*NyL + (i-1 - NyL) * (Nx+2)+j;
        end
        
        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        if (i == NyM+1) && (j < NxL + 2)
            col(ctr) = row(ctr) + (Nx+2);
        else 
            col(ctr) = row(ctr) + Nsml;
        end
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South 
        ctr = ctr + 1;
        row(ctr) = ROW;
        if (i == NyM+2) && (j > NxL+1)
            col(ctr) = row(ctr) - Nsml;
        else
            col(ctr) = row(ctr) - (Nx+2);
        end
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = -DE -DN -DW -DS + dx*dy/dt;
    
    end
end

% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+3 : Ny

    j_ctr = NxL+3; % Starting value for index on right of H
    
    for j = [2 : NxL+1, NxM+2 : Nx+1]
        dy = YU(i+1,j)-YU(i,j);
        yn = YV(i+1,j);
        yp = YV(i,j);
        ys = YV(i-1,j);
        dx = XU(i,j) - XU(i,j-1);
        xe = XV(i,j+1);
        xp = XV(i,j);
        xw = XV(i,j-1);
        
        if j < NxM+2
            ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j_ctr;
        end
        
        % East
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + 1;
        DE = - 1/Re * dy/(xe-xp)/2; 
        A(ctr) = DE;
        
        % North
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW + Nsml;
        DN = - 1/Re * dx/(yn-yp)/2;  
        A(ctr) = DN;

        % West
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - 1;
        DW = - 1/Re * dy/(xp-xw)/2;
        A(ctr) = DW;
        
        % South
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW - Nsml;
        DS = - 1/Re * dx/(yp-ys)/2;
        A(ctr) = DS;

        % Principal
        ctr = ctr + 1;
        row(ctr) = ROW;
        col(ctr) = ROW;
        A(ctr) = -DE -DN -DW -DS + dx*dy/dt;
    
    end
end

% Remove zero entries from row, col, A
A = nonzeros(A); 
row = nonzeros(row); 
col = nonzeros(col);

A = sparse(row, col, A, NNXY, NNXY);