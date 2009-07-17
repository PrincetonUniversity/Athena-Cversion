function [X,Y,status] = lineout_ijk(Grid,var,i,j,k)
% 
% lineout_ijk:  CREATE A 1D LINEOUT PLOT OF A GIVEN VARIABLE USING GRID
% COORDINATES.  EXACTLY ONE OF i, j, OR k MUST BE A VECTOR OF COORDINATES,
% THE REMAINING TWO BEING THE FIXED VALUES THROUGH WHICH TO TAKE THE
% LINEOUT.  
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/23/09

X = [];  Y = [];
status = 0;

[x,y,z] = ijk_to_xyz(Grid,i,j,k);

lenx = length(x);
leny = length(y);
lenz = length(z);

if (lenx>1 && leny==1 && lenz==1)
    X = x;
    Y = squeeze(var(:,j,k));
elseif (lenx==1 && leny>1 && lenz==1)
    X = y;
    Y = squeeze(var(i,:,k));
elseif (lenx==1 && leny==1 && lenz>1)
    X = z;
    Y = squeeze(var(i,j,:));
else
    status = -1;
    fprintf(2,'[lineout_ijk]:  Exactly one of i,j,k must be a vector!\n');
    return;
end;

return;