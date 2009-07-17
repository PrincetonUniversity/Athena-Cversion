global coordsys nx1 nx2 nx3 nvar nscalars ifgrav ndim dx1 dx2 dx3 ...
    x1min x1max x2min x2max x3min x3max gamma_1 iso_csound time dt ...
    base step path;

path = '/n/a2/askinner/nike/bin/cylblast/testing/';
% basename1 = 'CylBlast_B0';
basename1 = 'CylBlast_B0_Parallel_2_2';
basename2 = 'CylBlast_B0_Parallel_4_4';
step = 0;

minstep = 0;  maxstep = 2;
nstep = maxstep-minstep+1;

for step = minstep:maxstep
    filename1 = construct_filename(path,basename1,step);
    filename2 = construct_filename(path,basename2,step);
    [x1,x2,x3,d,M1,M2,M3,E,B1,B2,B3,status] = readbin(filename1);
    if (status==-1) 
        return;
    end;
    var_serial = M2;

    [x1,x2,x3,d,M1,M2,M3,E,B1,B2,B3,status] = readbin(filename2);
    if (status==-1) 
        return;
    end;
    var_parallel = M2;

    error = abs(var_serial-var_parallel);
    spy(error);
    max(max(error))

%     [X,Y,Z,status] = slice_xyz(error,x1,x2,0);
    [X,Y,Z,status] = slice_xyz(E,x1,x2,0);
    pseudocolor(Z);
    axis([1 1.5 x2min pi/4]);
    colorbar;

    pause;
end;