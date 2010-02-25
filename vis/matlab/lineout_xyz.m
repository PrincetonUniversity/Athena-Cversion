function [X,Y,status] = lineout_xyz(Grid,var,x,y,z,interp)
% 
% lineout_xyz:  CREATE A 1D LINEOUT PLOT OF A GIVEN VARIABLE USING SPATIAL
% COORDINATES.  EXACTLY ONE OF x, y, OR z MUST BE A VECTOR OF COORDINATES,
% THE REMAINING TWO BEING THE FIXED VALUES THROUGH WHICH TO TAKE THE
% LINEOUT.  THIS IS ESSENTIALLY A WRAPPER FOR lineout_ijk.
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  10/06/09

[i,j,k,onfacei,onfacej,onfacek] = xyz_to_ijk(Grid,x,y,z);

[X,Y,status] = lineout_ijk(Grid,var,i,j,k);

if (interp && onfacei)
    [X,Y2,status] = lineout_ijk(Grid,var,i+1,j,k);
    Y = 0.5*(Y+Y2);
end;
if (interp && onfacej)
    [X,Y2,status] = lineout_ijk(Grid,var,i,j+1,k);
    Y = 0.5*(Y+Y2);
end;
if (interp && onfacei)
    [X,Y2,status] = lineout_ijk(Grid,var,i,j,k+1);
    Y = 0.5*(Y+Y2);
end;

if (status == -1)
    printf(2,'[lineout_xyz]:  lineout_ijk returned in error!');
end;

return;