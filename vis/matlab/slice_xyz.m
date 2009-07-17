function [X,Y,Z,status] = slice_xyz(Grid,var,x,y,z)
% 
% slice_xyz:  CREATE A 2D SLICE OF A GIVEN VARIABLE USING SPATIAL 
% COORDINATES.  EXACTLY TWO OF x,y,z MUST BE VECTORS OF COORDINATES, THE 
% REMAINING ONE BEING THE FIXED VALUE THROUGH WHICH TO TAKE THE SLICE.
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/22/09

% TRANSFORM IN CYLINDRICAL CASE?
transform = true;

X = [];  Y = [];  Z = [];
status = 0;

% DETERMINE WHICH DIRECTION TO SLICE
lenx = length(x);  leny = length(y);  lenz = length(z);
if (lenx==1 && leny>1 && lenz>1) 
    dir = 1;
elseif (lenx>1 && leny==1 && lenz>1) 
    dir = 2;
elseif (lenx>1 && leny>1 && lenz==1) 
    dir = 3;
else
    status = -1;
    fprintf(2,'[slice_xyz]:  Exactly two of x,y,z must be vectors!\n');
    return;
end;

[i,j,k] = xyz_to_ijk(Grid,x,y,z);

% COMPUTE SLICE
switch (dir)
    case 1
        [X,Y] = meshgrid(y,z);
        Z = squeeze(var(i,:,:))';
    case 2
        [X,Y] = meshgrid(x,z);
        Z = squeeze(var(:,j,:))';
    case 3
        if ((Grid.coordsys == -2) && transform) % CYLINDRICAL
            [R,PHI] = meshgrid(x,y);
            [X,Y] = pol2cart(PHI,R);
        else
            [X,Y] = meshgrid(x,y);
        end;
        Z = squeeze(var(:,:,k))';
end;

return;