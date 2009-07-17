function map = my_maps(name,size)

upramp = linspace(0,1,size)';
downramp = linspace(1,0,size)';
hat = 1-abs(linspace(-1,1,size)');

switch (name)
    case 'xray'
        map = [downramp,downramp,downramp];
    case 'hot'
        map = [upramp,hat,downramp];
    otherwise
        fprintf(2,'[my_maps]:  %s is not a known colormap!\n',name);
        map = [];
end;

return;