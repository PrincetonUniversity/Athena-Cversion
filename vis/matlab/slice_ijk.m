function [X,Y,Z,status] = slice_ijk(Grid,var,i,j,k)
% 
% slice_ijk:  CREATE A 2D SLICE OF A GIVEN VARIABLE USING GRID COORDINATES.
% EXACTLY TWO OF i,j,k MUST BE VECTORS OF COORDINATES, THE REMAINING ONE 
% BEING THE FIXED VALUE THROUGH WHICH TO TAKE THE SLICE.  THIS IS
% ESSENTIALLY A WRAPPER FOR slice_xyz.
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/22/09

[x,y,z] = ijk_to_xyz(Grid,i,j,k);

[X,Y,Z,status = slice_xyz(Grid,var,x,y,z);

return;