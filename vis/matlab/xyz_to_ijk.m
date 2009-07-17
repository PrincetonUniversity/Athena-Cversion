function [i,j,k] = xyz_to_ijk(Grid,x,y,z)
% 
% xyz_to_ijk:  CONVERT A SPATIAL COORDINATE TO ITS CORRESPONDING GRID
% COORDINATE--SORT OF THE OPPOSITE OF ATHENA'S cc_pos.  
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/22/09

i = 1; j = 1; k = 1;
if (Grid.dx1>0) 
    i = round((x-Grid.x1min)./Grid.dx1-0.5)+1;
    i = max(i,1);  i = min(i,Grid.nx1);
end;
if (Grid.dx2>0)
    j = round((y-Grid.x2min)./Grid.dx2-0.5)+1;
    j = max(j,1);  j = min(j,Grid.nx2);
end;
if (Grid.dx3>0)
    k = round((z-Grid.x3min)./Grid.dx3-0.5)+1;
    k = max(k,1);  k = min(k,Grid.nx3);
end;