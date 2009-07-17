% 2D CASE
close all;

varname = 'd';
titlestr = 'DENSITY';
% varname = 'P';
% titlestr = 'PRESSURE';
% varname = 'Ekin';
% titlestr = 'SPECIFIC KINETIC ENERGY';
% varname = 'Emag';
% titlestr = 'MAGNETIC ENERGY';
% varname = 'A3';
% titlestr = 'MAGNETIC FIELD';
varname = 'normangmomflux'
titlestr = 'NORMALIZED ANGULAR MOMENTUM FLUX'

ncontours = 30;
maxstep = 30;

% CYLINDRICAL VERSION 
path = '/n/a2/askinner/nike/bin/cylrayleigh/';
% basename = 'CylRayleigh_1_95_Joined';
basename = 'CylRayleigh_1_99_Joined';
% basename = 'CylRayleigh_2_01_Joined';
% basename = 'CylRayleigh_2_05_Joined';

filename = construct_filename(path,basename,0);
[Grid,status] = init_grid(filename);

for step = 0:maxstep
    filename = construct_filename(path,basename,step);
    [time,dt,M1,status] = getvar(Grid,filename,'M1');
    [time,dt,M2,status] = getvar(Grid,filename,'M2');
    [time,dt,d,status] = getvar(Grid,filename,'d');
    [time,dt,P,status] = getvar(Grid,filename,'P');
    y = mean(M1.*M2./d,2)./mean(P,2);
    x = Grid.x1zones;
    y = y.*x/(0.5*(Grid.x1max+Grid.x1min));

%     y = 1 - log(mean(M2./d,2)/(2*pi))./log(x);
    plot(x,y);
%     axis([3 7 -.01 .01]);
    
%     % GET VARIABLE
%     [var,status] = getvar(Grid,Gas,varname);

    % PSEUDOCOLOR PLOT
%     plot1 = figure;
%     x = Grid.x1nodes;
%     y = Grid.x2nodes;
%     z = 0.0;
%     [X,Y,Z,status] = slice_xyz(Grid,var,x,y,z); 
%     my_pcolor(X,Y,Z);
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
    colormap(my_maps('hot',128));
    colorbar;
    titlestr = sprintf('step=%d, max=%f, min=%f',step,max(max(y)),min(min(y)));
    title(titlestr);

%     % CONTOUR PLOT
%     plot2 = figure;
%     x = Grid.x1zones;
%     y = Grid.x2zones;
%     z = 0.0;
%     [X,Y,Z,status] = slice_xyz(Grid,var,x,y,z); 
%     contour(X,Y,Z,ncontours,'k');
%     axis equal;
%     axis([1.5/sqrt(2)-.5 1.5/sqrt(2)+.5 1.5/sqrt(2)-.5 1.5/sqrt(2)+.5]);
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     title(titlestr);
% 
%     % LINEOUT PLOT
%     plot3 = figure;
%     x = Grid.x1zones;
%     y = pi/4.0;
%     z = 0.0;
%     [X,Y,status] = lineout_xyz(Grid,var,x,y,z);
%     plot(1:length(X)',Y,'ro');
%     title(titlestr);
    
    pause;
end;