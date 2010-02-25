function [X,Y,Z,status] = slice_xyz(Grid,var,x,y,z,interp)
% 
% slice_xyz:  CREATE A 2D SLICE OF A GIVEN VARIABLE USING SPATIAL 
% COORDINATES.  EXACTLY TWO OF x,y,z MUST BE VECTORS OF COORDINATES, THE 
% REMAINING ONE BEING THE FIXED VALUE THROUGH WHICH TO TAKE THE SLICE.  IF
% IT IS DETERMINED THAT THIS FIXED VALUE LIES BETWEEN GRID ZONES (E.G. FOR
% AN EVEN GRID), AN OPTIONAL INTERPOLATION FLAG MAY BE SET TO 1 AND THE
% AVERAGE OF THE TWO ADJACENT SLICES WILL BE RETURNED.  NOTE THAT IN THE
% CASE OF (ANTI-)SYMMETRY FOR AN EVEN GRID, THIS INTERPOLATION IS EXACT.
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  10/06/09

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

[i,j,k,onfacei,onfacej,onfacek] = xyz_to_ijk(Grid,x,y,z);

% COMPUTE SLICE
switch (dir)
    case 1
        [X,Y] = meshgrid(y,z);
        Z = squeeze(var(i,:,:))';
        if (onfacei && interp)
            Z = 0.5*(Z+squeeze(var(i+1,:,:))');
        end;
    case 2
        [X,Y] = meshgrid(x,z);
        Z = squeeze(var(:,j,:))';
        if (onfacej && interp)
            Z = 0.5*(Z+squeeze(var(:,j+1,:))');
        end;
    case 3
        if ((Grid.coordsys == -2) && transform) % CYLINDRICAL
            [R,PHI] = meshgrid(x,y);
            [X,Y] = pol2cart(PHI,R);
        else
            [X,Y] = meshgrid(x,y);
        end;
        Z = squeeze(var(:,:,k))';
        if (onfacek && interp)
            Z = 0.5*(Z+squeeze(var(:,:,k+1))');
        end;
end;

return;