function U = BoundCond(NumCells,GhostCells,U)
% This function assigns cell average values to the ghost cells 

% Input:
%  GhostCells: Index of ghost cells
%  NumCells: number of triangles in the initial mesh
%  U: cell average values on the mesh triangles

% Output:
%  U: Updated cell average values
%
IndT_r = GhostCells{1,1}; IndT_l = GhostCells{1,2}; IndT_u = GhostCells{1,3}; IndT_d = GhostCells{1,4};

% values on the left side
Len_r = length(IndT_r);
U(NumCells+1:NumCells+Len_r,1) = U(IndT_r,1);

% values on the right side
Len_l = length(IndT_l);
U(NumCells+Len_r+1:NumCells+Len_r+Len_l,1) = U(IndT_l,1);

% values on the up side 
Len_d = length(IndT_u);
U(NumCells+Len_r+Len_l+1:NumCells+Len_r+Len_l+Len_d,1) = U(IndT_u,1);

% values on the down side 
Len_u = length(IndT_d);
U(NumCells+Len_r+Len_l+Len_d+1:NumCells+Len_r+Len_l+Len_d+Len_u,1) = U(IndT_d,1);
