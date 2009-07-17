function map = my_maps(name,size)
% 
% my_maps:  CREATES A COLORMAP LIKE THOSE OF ATHENA THAT ARE
% (CONSPICUOUSLY) MISSING FROM THE MATLAB STANDARD COLORMAPS.
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/23/09

switch (name)
    case 'xray'
        map = flipud(gray(size));
    case 'hot'
        map = hsv2rgb([linspace(2/3,0,size)',ones(size,1),ones(size,1)]);\
    otherwise
        fprintf(2,'[my_maps]:  %s is not a known colormap!\n',name);
        map = jet(size);
end;

return;