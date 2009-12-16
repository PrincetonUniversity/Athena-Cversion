function [Rc Pc Zc T D Vr Vp Vz IntE Br Bp Bz] = ReadVTK(vtkfile, mhd)

% open file (big endian format)
fid = fopen(vtkfile,'r','b');

if( fid == -1 )
  return
end

%Read past header trash

line = fgetl(fid); % # vtk DataFile Version x.x
line = fgetl(fid); % comments
%Get Time @ Timestep
T = str2num(line(35:length(line)));
line = fgetl(fid); % BINARY
line = fgetl(fid); % DATASET STRUCTURED_POINTS

%Get dimensionality of data set
s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d').';
Nr = sz(1) - 1;
Np = sz(2) - 1;
Nz = sz(3) - 1;

%Get origin
line = fgetl(fid); % ORIGIN OX OY OZ
[Or] = sscanf(line, '%*s%e%e%e');
r0 = Or(1);
p0 = Or(2);
z0 = Or(3);

%Get spacing
line = fgetl(fid); % SPACING SX SY SZ
[Dx] = sscanf(line, '%*s%e%e%e');
dr = Dx(1);
dp = Dx(2);
dz = Dx(3);

%Calculate cell centered spatial values
Ri = r0 + (dr)*(0:Nr); Rc = Ri + (dr/2); Rc = Rc(1:Nr);
Pi = p0 + (dp)*(0:Np); Pc = Pi + (dp/2); Pc = Pc(1:Np);
Zi = z0 + (dz)*(0:Nz); Zc = Zi + (dz/2); Zc = Zc(1:Nz);

line = fgetl(fid); % POINT_DATA NXNYNZ

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
line = fgetl(fid);

Elts = Nr*Np*Nz; Dat = 'float';
%Read density
    D = fread(fid,Elts,Dat);
    D = reshape(D,Nr,Np,Nz);
    
s = fgetl(fid); 
s = fgetl(fid); 
   
%Read Mr, convert to Vr
    Vr = fread(fid,Elts,Dat);
    Vr = reshape(Vr,Nr,Np,Nz);
    Vr = Vr./D;
    
s = fgetl(fid); 
s = fgetl(fid); 

%Read Mphi, convert to Vp
    Vp = fread(fid,Elts,Dat);
    Vp = reshape(Vp,Nr,Np,Nz);
    Vp = Vp./D;
    
s = fgetl(fid); 
s = fgetl(fid); 

%Read Mz, convert to Vz
    Vz = fread(fid,Elts,Dat);
    Vz = reshape(Vz,Nr,Np,Nz);
    Vz = Vz./D;
    
s = fgetl(fid); %disp(s)
s = fgetl(fid); %disp(s)

%Read IntE
    IntE = fread(fid,Elts,Dat);
    IntE = reshape(IntE,Nr,Np,Nz);
    
%Read magnetic variables if necessary (mhd == 1)
if (mhd == 1 )
    s = fgetl(fid); %disp(s)
    s = fgetl(fid); %disp(s)
    %Read Br
    Br = fread(fid,Elts,Dat);
    Br = reshape(Br,Nr,Np,Nz);
    
    s = fgetl(fid); %disp(s)
    s = fgetl(fid); %disp(s)
    %Read Bp
    Bp = fread(fid,Elts,Dat);
    Bp = reshape(Bp,Nr,Np,Nz);
    
    s = fgetl(fid); %disp(s)
    s = fgetl(fid); %disp(s)
    %Read Bz
    Bz = fread(fid,Elts,Dat);
    Bz = reshape(Bz,Nr,Np,Nz);
    
else
    Br = D;
    Bp = D;
    Bz = D;
end

fclose(fid);
return;

