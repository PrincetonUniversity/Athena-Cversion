function [Gas,status] = init_gas(Grid,first,last)

status = 0;

n = 1;
nx1 = Grid.nx1;
nx2 = Grid.nx2;
nx3 = Grid.nx3;

for step = first:last
    Gas(n).step = step;
    Gas(n).time = 0;
    Gas(n).dt = 0;
    Gas(n).d  = squeeze(zeros(nx1,nx2,nx3));
    Gas(n).M1 = squeeze(zeros(nx1,nx2,nx3));
    Gas(n).M2 = squeeze(zeros(nx1,nx2,nx3));
    Gas(n).M3 = squeeze(zeros(nx1,nx2,nx3));
    if (Grid.adiabatic)
        Gas(n).E = squeeze(zeros(nx1,nx2,nx3));
    end;
    if (Grid.mhd)
        Gas(n).B1 = squeeze(zeros(nx1,nx2,nx3));
        Gas(n).B2 = squeeze(zeros(nx1,nx2,nx3));
        Gas(n).B3 = squeeze(zeros(nx1,nx2,nx3));
    end;
    n = n + 1;
end;

return;