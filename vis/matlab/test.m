global nx1 nx2 nx3 nvar nscalars ifgrav ndim dx1 dx2 dx3 ...
    x1min x1max x2min x2max x3min x3max gamma_1 iso_csound time dt ...
    base step path;

path = '/n/a2/askinner/nike/bin/cylwindrotb/2D/allregions/';
basename = 'CylWindRotB_128';
step = 0;

minstep = 0;  maxstep = 50;
nstep = maxstep-minstep+1;

n = 0;
data_d = zeros(nx1,nx2,nx3,nstep);

for step = minstep:maxstep
    filename = construct_filename(path,basename,step);
    [x1,x2,x3,d,M1,M2,M3,E,B1,B2,B3,status] = readbin(filename);
    if (status==-1) 
        return;
    end;
    n = n + 1;
    data(:,:,:,n) = M1;
    
    error = data(:,:,:,n)-data(:,:,:,1);
    [X,Y,status] = lineout_xyz(error,x1,0,0);

    linespec = '-o';
    xlbl = 'R';
    ylbl = 'Error_d';
    
    plot(X,Y,linespec);
    titlestr = sprintf('%s\t\tStep %d,\tTime %f',base,step,time);
    titlestr = strrep(titlestr,'_','\_');
    xlbl = strrep(xlbl,'_','\_');
    ylbl = strrep(ylbl,'_','\_');
    title(titlestr);
    xlabel(xlbl);
    ylabel(ylbl);
    axis tight;
    
    pause;
end;