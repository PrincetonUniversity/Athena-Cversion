function [x,y,z] = ijk_to_xyz(Grid,i,j,k)
% 
% ijk_to_xyz:  CONVERT A GRID COORDINATE TO ITS CORRESPONDING SPATIAL
% (NOT NECESSARILY GEOMETRIC) CENTER COORDINATE--SORT OF LIKE ATHENA'S 
% cc_pos.  
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/22/09

x = Grid.x1min + (i-1+0.5)*Grid.dx1;
y = Grid.x2min + (j-1+0.5)*Grid.dx2;
z = Grid.x3min + (k-1+0.5)*Grid.dx3;