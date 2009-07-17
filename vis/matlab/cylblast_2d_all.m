% 2D CASE
close all;

% varname = 'P';
% titlestr = 'PRESSURE';
% varname = 'Ekin';
% titlestr = 'SPECIFIC KINETIC ENERGY';
% varname = 'Emag';
% titlestr = 'MAGNETIC ENERGY';
% varname = 'A3';
% titlestr = 'MAGNETIC FIELD';

ncontours = 30;
scrsz = get(0,'ScreenSize');
figwidth = scrsz(3)/2;
figheight = scrsz(4);
scl = .04;

% CYLINDRICAL VERSION 
path = '/n/a2/askinner/nike/bin/cylblast/2D/';
basename = 'CylBlast_B0_Joined';
step = 4;

filename = construct_filename(path,basename,0);
[Grid_cyl,status] = init_grid(filename);
filename = construct_filename(path,basename,step);


varname = 'd';
titlestr = 'DENSITY';
% GET VARIABLE
[time,dt,var,status] = getvar(Grid_cyl,filename,varname);

% PSEUDOCOLOR PLOT
plot1 = figure('Position',[1 1 figwidth figheight]);
h = subplot(3,2,1);
p = get(h,'pos');
p(1) = p(1) - scl/2;
p(2) = p(2) - scl/2;
p(3) = p(3) + scl;
p(4) = p(4) + scl;
set(h,'pos',p);
x = Grid_cyl.x1nodes;
y = Grid_cyl.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
my_pcolor(X,Y,Z);
axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(my_maps('hot',128));
colorbar;
title(titlestr);

% % CONTOUR PLOT
% plot2 = figure('Position',[1 1 figwidth figheight]);
% subplot(3,2,1,'align');
% x = Grid_cyl.x1zones;
% y = Grid_cyl.x2zones;
% z = 0.0;
% [X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
% contour(X,Y,Z,ncontours,'k');
% axis equal;
% axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% title(titlestr);
% 
% % LINEOUT PLOT
% plot3 = figure('Position',[1 1 figwidth figheight]);
% subplot(3,1,1,'align');
% x = Grid_cyl.x1zones;
% y = pi/4.0;
% z = 0.0;
% [X,Y,status] = lineout_xyz(Grid_cyl,var,x,y,z);
% plot(1:length(X)',Y,'ro');
% title(titlestr);



varname = 'P';
titlestr = 'PRESSURE';
% GET VARIABLE
[time,dt,var,status] = getvar(Grid_cyl,filename,varname);

% PSEUDOCOLOR PLOT
figure(plot1);
h = subplot(3,2,3);
p = get(h,'pos');
p(1) = p(1) - scl/2;
p(2) = p(2) - scl/2;
p(3) = p(3) + scl;
p(4) = p(4) + scl;
set(h,'pos',p);
x = Grid_cyl.x1nodes;
y = Grid_cyl.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
my_pcolor(X,Y,Z);
axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(my_maps('hot',128));
colorbar;
title(titlestr);

% % CONTOUR PLOT
% figure(plot2);
% subplot(3,2,3,'align');
% x = Grid_cyl.x1zones;
% y = Grid_cyl.x2zones;
% z = 0.0;
% [X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
% contour(X,Y,Z,ncontours,'k');
% axis equal;
% axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% title(titlestr);
% 
% % LINEOUT PLOT
% figure(plot3);
% subplot(3,1,2,'align');
% x = Grid_cyl.x1zones;
% y = pi/4.0;
% z = 0.0;
% [X,Y,status] = lineout_xyz(Grid_cyl,var,x,y,z);
% plot(1:length(X)',Y,'ro');
% title(titlestr);


varname = 'Ekin';
titlestr = 'SPECIFIC KINETIC ENERGY';
% GET VARIABLE
[time,dt,var,status] = getvar(Grid_cyl,filename,varname);

% PSEUDOCOLOR PLOT
figure(plot1);
h = subplot(3,2,5);
p = get(h,'pos');
p(1) = p(1) - scl/2;
p(2) = p(2) - scl/2;
p(3) = p(3) + scl;
p(4) = p(4) + scl;
set(h,'pos',p);
x = Grid_cyl.x1nodes;
y = Grid_cyl.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
my_pcolor(X,Y,Z);
axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(my_maps('hot',128));
colorbar;
title(titlestr);

% % CONTOUR PLOT
% figure(plot2);
% subplot(3,2,5,'align');
% x = Grid_cyl.x1zones;
% y = Grid_cyl.x2zones;
% z = 0.0;
% [X,Y,Z,status] = slice_xyz(Grid_cyl,var,x,y,z); 
% contour(X,Y,Z,ncontours,'k');
% axis equal;
% axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% title(titlestr);
% 
% % LINEOUT PLOT
% figure(plot3);
% subplot(3,1,3,'align');
% x = Grid_cyl.x1zones;
% y = pi/4.0;
% z = 0.0;
% [X,Y,status] = lineout_xyz(Grid_cyl,var,x,y,z);
% plot(1:length(X)',Y,'ro');
% title(titlestr);
% 



% CARTESIAN VERSION
path = '/n/a2/askinner/nike/bin/cylblast/2D/';
basename = 'Blast_B0_Joined';
step = 4;

filename = construct_filename(path,basename,0);
[Grid_cart,status] = init_grid(filename);
filename = construct_filename(path,basename,step);


varname = 'd';
titlestr = 'DENSITY';
% GET VARIABLE
[time,dt,var,status] = getvar(Grid_cart,filename,varname);

% PSEUDOCOLOR PLOT
figure(plot1);
h = subplot(3,2,2);
p = get(h,'pos');
p(1) = p(1) - scl/2;
p(2) = p(2) - scl/2;
p(3) = p(3) + scl;
p(4) = p(4) + scl;
set(h,'pos',p);
x = Grid_cart.x1nodes;
y = Grid_cart.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
my_pcolor(X,Y,Z);
axis([-.5 .5 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(my_maps('hot',128));
colorbar;
title(titlestr);

% % CONTOUR PLOT
% figure(plot2);
% subplot(3,2,2,'align');
% x = Grid_cart.x1zones;
% y = Grid_cart.x2zones;
% z = 0.0;
% [X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
% contour(X,Y,Z,ncontours,'k');
% axis equal;
% axis([-.5 .5 -.5 .5]);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% title(titlestr);
% 
% % LINEOUT PLOT
% figure(plot3);
% subplot(3,1,1,'align');
% x = Grid_cart.x1zones;
% y = 0.0;
% z = 0.0;
% [X,Y,status] = lineout_xyz(Grid_cart,var,x,y,z);
% hold on;
% plot(1:length(X)',Y);
% axis tight;
% hold off;


varname = 'P';
titlestr = 'PRESSURE';
% GET VARIABLE
[time,dt,var,status] = getvar(Grid_cart,filename,varname);

% PSEUDOCOLOR PLOT
figure(plot1);
h = subplot(3,2,4);
p = get(h,'pos');
p(1) = p(1) - scl/2;
p(2) = p(2) - scl/2;
p(3) = p(3) + scl;
p(4) = p(4) + scl;
set(h,'pos',p);
x = Grid_cart.x1nodes;
y = Grid_cart.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
my_pcolor(X,Y,Z);
axis([-.5 .5 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(my_maps('hot',128));
colorbar;
title(titlestr);

% % CONTOUR PLOT
% figure(plot2);
% subplot(3,2,4,'align');
% x = Grid_cart.x1zones;
% y = Grid_cart.x2zones;
% z = 0.0;
% [X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
% contour(X,Y,Z,ncontours,'k');
% axis equal;
% axis([-.5 .5 -.5 .5]);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% title(titlestr);
% 
% % LINEOUT PLOT
% figure(plot3);
% subplot(3,1,2,'align');
% x = Grid_cart.x1zones;
% y = 0.0;
% z = 0.0;
% [X,Y,status] = lineout_xyz(Grid_cart,var,x,y,z);
% hold on;
% plot(1:length(X)',Y);
% axis tight;
% hold off;


varname = 'Ekin';
titlestr = 'SPECIFIC KINETIC ENERGY';
% GET VARIABLE
[time,dt,var,status] = getvar(Grid_cart,filename,varname);

% PSEUDOCOLOR PLOT
figure(plot1);
h = subplot(3,2,6);
p = get(h,'pos');
p(1) = p(1) - scl/2;
p(2) = p(2) - scl/2;
p(3) = p(3) + scl;
p(4) = p(4) + scl;
set(h,'pos',p);
x = Grid_cart.x1nodes;
y = Grid_cart.x2nodes;
z = 0.0;
[X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
my_pcolor(X,Y,Z);
axis([-.5 .5 -.5 .5]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
colormap(my_maps('hot',128));
colorbar;
title(titlestr);

% % CONTOUR PLOT
% figure(plot2);
% subplot(3,2,6,'align');
% x = Grid_cart.x1zones;
% y = Grid_cart.x2zones;
% z = 0.0;
% [X,Y,Z,status] = slice_xyz(Grid_cart,var,x,y,z); 
% contour(X,Y,Z,ncontours,'k');
% axis equal;
% axis([-.5 .5 -.5 .5]);
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% title(titlestr);
% 
% % LINEOUT PLOT
% figure(plot3);
% subplot(3,1,3,'align');
% x = Grid_cart.x1zones;
% y = 0.0;
% z = 0.0;
% [X,Y,status] = lineout_xyz(Grid_cart,var,x,y,z);
% hold on;
% plot(1:length(X)',Y);
% axis tight;
% hold off;