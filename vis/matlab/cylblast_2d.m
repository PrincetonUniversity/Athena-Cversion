% 2D CASE
close all;

varname = 'd';
titlestr = 'DENSITY';

% varname = 'P';
% titlestr = 'PRESSURE';

% varname = 'SpecEkin';
% titlestr = 'SPECIFIC KINETIC ENERGY';

% varname = 'Emag';
% titlestr = 'MAGNETIC ENERGY';

% varname = 'A3';
% titlestr = 'MAGNETIC FIELD';

ncontours = 30;

% CYLINDRICAL VERSION 
% path = '/n/a2/askinner/nike/bin/cylblast/2D/';
% basename = 'CylBlast_B10_Joined';
path = '/n/a2/askinner/svnathena/';
basename = 'CylBlast_B10';
step = 4;

filename = construct_filename(path,basename,step);
[Grid_cyl,status] = init_grid(filename);

% GET VARIABLE
filename = construct_filename(path,basename,step);
[time,dt,var,status] = getvar(Grid_cyl,filename,varname);

% PSEUDOCOLOR PLOT
plot1 = figure;
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 4 4]);
x = Grid_cyl.x1nodes;
y = Grid_cyl.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
my_pcolor(X,Y,Z);
% axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
axis([1 2 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(mymaps('hot',128));
% colorbar;
title(titlestr);

% CONTOUR PLOT
plot2 = figure;
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 4 4]);
x = Grid_cyl.x1zones;
y = Grid_cyl.x2zones;
z = 0.5;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
contour(X,Y,Z,ncontours,'k');
axis equal;
% axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
axis([1 2 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title(titlestr);

% LINEOUT PLOT
plot3 = figure;
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 10 4]);
x = Grid_cyl.x1zones;
y = 0.0;
z = 0.0;
[X,Y,status] = lineout_xyz(Grid_cyl,var,x,y,z);
plot(1:length(X)',Y,'ro');
title(titlestr);


% ncontours = 50;


% CARTESIAN VERSION 
path = '/n/a2/askinner/svnathena/';
basename = 'Blast_B10';
step = 4;

filename = construct_filename(path,basename,step);
[Grid_cart,status] = init_grid(filename);

% GET VARIABLE
filename = construct_filename(path,basename,step);
[time,dt,var,status] = getvar(Grid_cart,filename,varname);

% PSEUDOCOLOR PLOT
plot4 = figure;
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 4 4]);
x = Grid_cart.x1nodes;
y = Grid_cart.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
my_pcolor(X,Y,Z);
colormap(mymaps('hot',128));
axis([-.5 .5 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
% colorbar;
title(titlestr);

% CONTOUR PLOT
plot5 = figure;
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 4 4]);
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

% LINEOUT PLOT
figure(plot3);
x = Grid_cart.x1zones;
y = 0.0;
z = 0.0;
[X,Y,status] = lineout_xyz(Grid_cart,var,x,y,z);
hold on;
plot(1:length(X)',Y);
axis tight;
hold off;