
%ATH_READBIN    Read in Athena .bin files
% 
%   [TIME,DT,VAR,STATUS] = ATH_READBIN(GRID,FILENAME,VARNAME) reads .bin
%   file at location FILENAME and returns the desired variable indicated by
%   the string VARNAME.  If the flag consistency_check is set to true, then
%   a check is performed to ensure that the same grid metadata is used in
%   each opened file.
%
%   The only supported strings for VARNAME are
%       'd'             - density
%       'M1','M2','M3'  - momentum density components
%       'E'             - total energy density
%       'B1','B2','B3'  - magnetic field components
%       'Phi'           - gravitational potential
%       'Er'            - radiation energy density
%       'F1','F2','F3'  - radiation flux components
%
%   For a longer list of derived variables, see ATH_GETVAR
%
%   AUTHOR:  Aaron Skinner
%   LAST MODIFIED:  11/17/2010
function [time,dt,var,status] = ath_readbin(Grid,filename,varname)

% realtype = 'single';
realtype = 'double';
consistency_check = false;
time = 0.0;
dt = 0.0;
var = [];
status = 0;

% PARSE FILENAME AND TEST TO SEE IF .bin
[path,basename,step,ext] = ath_parse_filename(filename);
if (~strcmp(ext,'.bin'))
    fprintf(2,'[ath_readbin]:  %s is not a .bin file!\n', filename);
    status = -1;
    return;
end;

% OPEN FILE FOR READING
[fid, message] = fopen(filename,'r');
if (fid==-1)
    fprintf(2,'[ath_readbin]:  %s could not be opened!\n', filename);
    fprintf(2,'%s', message);
    status = -1;
    return;
end;

if (consistency_check)
    % READ COORDINATE SYSTEM INFORMATION
    coordsys = fread(fid,1,'int');

    % READ NUMBER OF DATA POINTS
    dat = fread(fid,7,'int');
    nx1      = dat(1);
    nx2      = dat(2);
    nx3      = dat(3);
    nvar     = dat(4);
    nscalars = dat(5);
    selfg    = dat(6);
    ifpart   = dat(7);

    % READ (Gamma-1), ISOTHERMAL SOUND SPEED, TIME, AND dt
    dat = fread(fid,4,realtype);
    gamma_1    = dat(1);
    iso_csound = dat(2);
    time       = dat(3);
    dt         = dat(4);

    % READ X1,X2,X3 COORDINATES
    x1 = fread(fid,nx1,realtype);
    x2 = fread(fid,nx2,realtype);
    x3 = fread(fid,nx3,realtype);

    % CHECK FOR CONSISTENCY
    warn = false;
    warn = warn && (coordsys   == Grid.coordsys  );
    warn = warn && (nx1        == Grid.nx1       );
    warn = warn && (nx2        == Grid.nx2       );
    warn = warn && (nx3        == Grid.nx3       );
    warn = warn && (nvar       == Grid.nvar      );
    warn = warn && (nscalars   == Grid.nscalars  );
    warn = warn && (ifselfg    == Grid.selfg     );
    warn = warn && (ifpart     == Grid.particles );
    warn = warn && (gamma_1    == Grid.gamma_1   );
    warn = warn && (iso_csound == Grid.iso_csound);
    if (warn)
        fprintf(2, ...
            '[ath_readbin]:  %s failed consistency check!\n',filename);
        status = -1;
        return;
    else
        fprintf('CONSISTENCY CHECK PASSED!\n');
    end;
end;

% COMPUTE OFFSET FOR DESIRED VARIABLE
nx1 = Grid.nx1;
nx2 = Grid.nx2;
nx3 = Grid.nx3;
N = nx1*nx2*nx3;

eskip = (Grid.adiabatic > 0);  % SKIP ENERGY IF ADIABATIC DEFINED
bskip = 3*(Grid.mhd > 0);  % SKIP B-FIELDS IF MHD DEFINED
gskip = (Grid.selfg > 0);  % SKIP PHI IF SELFG DEFINED
switch(varname)
    case 'd'
        skip = 0;  % SKIP NONE, d IS FIRST
    case 'M1'
        skip = 1;  % SKIP 1 FOR d
    case 'M2'
        skip = 2;  % SKIP 2 FOR d,M1
    case 'M3'
        skip = 3;  % SKIP 3 FOR d,M1,M2
    case 'E'
        if (Grid.adiabatic)
            skip = 4;  % SKIP 4 FOR d,M1,M2,M3
        end;
    case 'B1'
        if (Grid.mhd)
            skip = eskip + 4;  % SKIP 4 FOR d,M1,M2,M3
        end;
    case 'B2'
        if (Grid.mhd)
            skip = eskip + 5;  % SKIP 5 FOR d,M1,M2,M3,B1
        end;
    case 'B3'
        if (Grid.mhd)
            skip = eskip + 6;  % SKIP 6 FOR d,M1,M2,M3,B1,B2
        end;
    case 'Phi'
        if (Grid.selfg)
            skip = eskip + bskip + 4;
        end;
    case 'Er'
        if (Grid.rad)
            skip = eskip + bskip + gskip + 4;
        end;
    case 'F1'
        if (Grid.rad)
            skip = eskip + bskip + gskip + 5;
        end;
    case 'F2'
        if (Grid.rad)
            skip = eskip + bskip + gskip + 6;
        end;
    case 'F3'
        if (Grid.rad)
            skip = eskip + bskip + gskip + 7;
        end;
    otherwise
        status = -1;
        fprintf('[ath_readbin]:  %s is not a valid variable!\n',varname);
        return;
end;

% READ IN time, dt
fseek(fid,Grid.time_offset,'bof');
time = fread(fid,1,realtype);
dt   = fread(fid,1,realtype);

% SET THE FILE POINTER TO THE BEGINNING OF THE DATA
offset = Grid.data_offset + skip*N*ath_sizeof(realtype);
status = fseek(fid,offset,'bof');
if (status == -1)
    fprintf('[ath_readbin]:  %s does not contain the requested data!\n',...
        filename);
    return;
end;

% READ IN CELL-CENTERED DATA
% Note:  The dimensions are permuted so that they match MATLAB's meshgrid
var = permute(reshape(fread(fid,N,realtype),nx1,nx2,nx3),[2 1 3]);

% CLOSE FILE
status = fclose(fid);
if (status == -1)
    fprintf('[ath_readbin]:  %s could not be closed!\n', filename);
end;

return;