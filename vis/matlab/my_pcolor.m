function status = my_pcolor(X,Y,C);
% 
% pseudocolor:  CREATE A PSEUDOCOLOR PLOT OF THE 2D DATA CONTAINED IN C.
% CONVERT COORDINATES IF NECESSARY. 
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/23/09

status = 0;

% CHECK COMPATIBILITY OF ARGUMENTS
[nx1,nx2] = size(X);
[ny1,ny2] = size(Y);
[nc1,nc2] = size(C);
if ~(nx1==ny1 && nx2==ny2)
    status = -1;
    fprintf(2,'[my_pcolor]:  X and Y must be the same size!\n');
    return;
end;
if ~((nc1==nx1 || nc1==nx1-1) && (nc2==nx2 || nc2==nx2-1))
    status = -1;
    fprintf(2,'[my_pcolor]:  C is not of compatible size!\n');
    return;
end;
    
% MAKE PSEUDOCOLOR PLOT
surf(X,Y,zeros(nx1,nx2),C);
view(2);
% axis equal;
% axis tight;
shading flat;

return;