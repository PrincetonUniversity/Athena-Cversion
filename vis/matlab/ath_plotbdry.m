%ATH_PLOTBDRY    Plot the 2D boundary of given grid.
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  2/1/2010
function ath_plotbdry(grid)

X = grid.x1nodes;
Y = grid.x2nodes;
Xbdry = [X;MAX(X)*ones(size(Y));flipud(X);MIN(X)*ones(size(Y))];
Ybdry = [MIN(Y)*ones(size(X));Y;MAX(Y)*ones(size(X));flipud(Y)];

if (grid.coordsys == -2)
    [Xbdry,Ybdry] = pol2cart(Ybdry,Xbdry);
end;

hold on;
plot(Xbdry,Ybdry,'k');
hold off;