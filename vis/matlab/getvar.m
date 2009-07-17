function [time,dt,var,status] = getvar(Grid,filename,varname)
% 
% getvar:  TAKES IN A GRID DATA STRUCTURE AND STRING ARGUMENTS FOR THE
% FILENAME AND DESIRED VARIABLE AND RETURNS THAT VARIABLE, IF IT EXISTS.  
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/29/09

status = 0;
time = 0.0;
dt = 0.0;
var = [];

switch(varname)
    case {'d','M1','M2','M3','E','B1','B2','B3'}
        [time,dt,var,status] = readbin(Grid,filename,varname);
        return;
    case {'Eint','Ekin','SpecEkin','Emag','P'}
        [time,dt,var,status] = readbin(Grid,filename,'M1');
        Ekin = var.^2;
        [time,dt,var,status] = readbin(Grid,filename,'M2');
        Ekin = Ekin + var.^2;
        [time,dt,var,status] = readbin(Grid,filename,'M3');
        Ekin = 0.5*(Ekin + var.^2);
        [time,dt,var,status] = readbin(Grid,filename,'d');
        Ekin = Ekin./var;
        if (strcmp(varname,'Ekin'))
            var = Ekin;
            return;
        end;
        if (strcmp(varname,'SpecEkin'))
            var = Ekin./var;
            return;
        end;
        Eint = 0.0;
        Emag = 0.0;
        if (Grid.mhd)
            [time,dt,var,status] = readbin(Grid,filename,'B1');
            Emag = var.^2;
            [time,dt,var,status] = readbin(Grid,filename,'B2');
            Emag = Emag + var.^2;
            [time,dt,var,status] = readbin(Grid,filename,'B3');
            Emag = 0.5*(Emag + var.^2);
            if (varname=='Emag')
                var = Emag;
                return;
            end;
        end;
        if (Grid.adiabatic)
            [time,dt,var,status] = readbin(Grid,filename,'E');
            Eint = var - Ekin - Emag;
            if (varname=='Eint')
                var = Eint;
                return;
            elseif (varname=='P')
                var = Eint*Grid.gamma_1;
                return;
            end;
        end;
    case 'A2'
        if (Grid.mhd)
            [time,dt,B1,status] = readbin(Grid,filename,'B1');
            [time,dt,B3,status] = readbin(Grid,filename,'B3');
            nx1 = Grid.nx1;
            nx2 = Grid.nx2;
            nx3 = Grid.nx3;
            dx1 = Grid.dx1;
            dx3 = Grid.dx3;
            x1 = Grid.x1zones;
            var = zeros(nx1,nx2,nx3);
            for j = 1:nx2
                for i = 2:nx1
                    x1i = 0.5*(x1(i) + x1(i-1));
                    B3i = (x1(i)*B3(i,j,1) + x1(i-1)*B3(i-1,j,1))/(2*x1i);
                    var(i,j,1) = (x1(i-1)*var(i-1,j,1)+x1i*dx1*B3i)/x1(i);
                end;
                for k = 2:nx3
                    for i = 1:nx1
                        B1i = 0.5*(B1(i,j,k) + B1(i,j,k-1));
                        var(i,j,k) = var(i,j,k-1) - B1i*dx3;
                    end;
                end;
            end;
            return;
        end;
    case 'A3'
        if (Grid.mhd)
            [time,dt,B1,status] = readbin(Grid,filename,'B1');
            [time,dt,B2,status] = readbin(Grid,filename,'B2');
            nx1 = Grid.nx1;
            nx2 = Grid.nx2;
            nx3 = Grid.nx3;
            dx1 = Grid.dx1;
            dx2 = Grid.dx2;
            x1 = Grid.x1zones;
            var = zeros(nx1,nx2,nx3);
            for k = 1:nx3

%                 % INTEGRATE X1, THEN X2
%                 for i = 2:nx1
%                     B2i = 0.5*(B2(i-1,1,k) + B2(i,1,k));
%                     var(i,1,k) = var(i-1,1,k) - B2i*dx1;
%                 end;
%                 for j = 2:nx2
%                     for i = 1:nx1
%                         B1i = 0.5*(B1(i,j-1,k) + B1(i,j,k));
%                         if (Grid.coordsys==-1)  % CARTESIAN
%                             var(i,j,k) = var(i,j-1,k) + B1i*dx2;
%                         elseif (Grid.coordsys==-2)  % CYLINDRICAL
%                             var(i,j,k) = var(i,j-1,k) + B1i*x1(i)*dx2;
%                         end;
%                     end;
%                 end;

                % INTEGRATE X2, THEN X1
                for j = 2:nx2
                    B1i = 0.5*(B1(1,j-1,k) + B1(1,j,k));
                    if (Grid.coordsys==-1)  % CARTESIAN
                        var(1,j,k) = var(1,j-1,k) + B1i*dx2;
                    elseif (Grid.coordsys==-2)  % CYLINDRICAL
                        var(1,j,k) = var(1,j-1,k) + B1i*x1(1)*dx2;
                    end;
                end;
                for i = 2:nx1
                    for j = 1:nx2
                        B2i = 0.5*(B2(i-1,j,k) + B2(i,j,k));
                        var(i,j,k) = var(i-1,j,k) - B2i*dx1;
                    end;
                end;
                
                
            end;
            
            return;
        end;
end;

status = -1;
fprintf('[getvar]:  %s is not a valid variable!\n',varname);
return;