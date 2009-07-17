% 3D CASE
close all;

% varname = 'd';
% titlestr = 'DENSITY';

varname = 'P';
titlestr = 'PRESSURE';

% varname = 'Ekin';
% titlestr = 'SPECIFIC KINETIC ENERGY';

% varname = 'Emag';
% titlestr = 'MAGNETIC ENERGY';

% varname = 'A3';
% titlestr = 'MAGNETIC FIELD';

ncontours = 30;

% CYLINDRICAL VERSION 
path = '/n/a2/askinner/nike/bin/cylblast/3D/';
basename = 'CylBlast_B0_Joined';
step = 20;

filename = construct_filename(path,basename,20);
[Grid_cyl,status] = init_grid(filename);

% GET VARIABLE
filename = construct_filename(path,basename,step);
[time,dt,var,status] = getvar(Grid_cyl,filename,varname);

% PSEUDOCOLOR PLOT
plot1 = figure;
x = Grid_cyl.x1nodes;
y = Grid_cyl.x2nodes;
z = 0.5;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
my_pcolor(X,Y,Z);
axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(my_maps('hot',128));
colorbar;
title(titlestr);

% CONTOUR PLOT
plot2 = figure;
x = Grid_cyl.x1zones;
y = Grid_cyl.x2zones;
z = 0.5;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
contour(X,Y,Z,ncontours,'k');
axis equal;
axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title(titlestr);



% CARTESIAN VERSION 
path = '/n/a2/askinner/nike/bin/cylblast/3D/';
basename = 'Blast_B0_Joined';
step = 20;

filename = construct_filename(path,basename,20);
[Grid_cart,status] = init_grid(filename);

% GET VARIABLE
filename = construct_filename(path,basename,step);
[time,dt,var,status] = getvar(Grid_cart,filename,varname);

% PSEUDOCOLOR PLOT
plot3 = figure;
x = Grid_cart.x1nodes;
y = Grid_cart.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
my_pcolor(X,Y,Z);
colormap(my_maps('hot',128));
axis([-.5 .5 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colorbar;
title(titlestr);

% CONTOUR PLOT
plot4 = figure;
x = Grid_cart.x1zones;
y = Grid_cart.x2zones;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
contour(X,Y,Z,ncontours,'k');
axis equal;
axis([-.5 .5 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title(titlestr);