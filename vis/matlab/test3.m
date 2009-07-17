% 2D CASE
close all;

varname = 'A3';
titlestr = 'A3';

% CYLINDRICAL VERSION 
path = '/n/a2/askinner/nike/bin/cylblast/2D/';
basename = 'CylBlast_B1';
step = 5;

filename = construct_filename(path,basename,0);
[Grid_cyl,status] = init_grid(filename);

filename = construct_filename(path,basename,step);
[Gas_cyl, status] = readbin(Grid_cyl,filename);

% GET VARIABLE
[var,status] = getvar(Grid_cyl,Gas_cyl,varname);

% CONTOUR PLOT
plot2 = figure;
x = Grid_cyl.x1zones;
y = Grid_cyl.x2zones;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
contour(X,Y,Z,100,'k');
axis square;
axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
xlabel('x');
ylabel('y');
title(titlestr);
