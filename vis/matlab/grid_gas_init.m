function [Grid,Gas,status] = grid_gas_init(first,last)

n = 1;

for step = first:last
    Grid(n).step = step;
    Grid(n).coordsys = 0;
    Grid(n).nx1 = 0;
    Grid(n).nx2 = 0;
    Grid(n).nx3 = 0;
    
    