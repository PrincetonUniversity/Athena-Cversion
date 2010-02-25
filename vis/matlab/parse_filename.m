function [path,basename,step,ext] = parse_filename(filename)
% 
% parse_filename:  BREAK UP A FULL-PATH FILENAME INTO ITS COMPONENT PARTS
% TO CHECK THE EXTENSION, MAKE IT MORE READABLE, AND EXTRACT THE STEP
% NUMBER.
% 
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/22/09

[path, file, ext, versn] = fileparts(filename);
path = strcat(path,'/');
[basename,remain] = strtok(file,'.');
while 1
    [token,remain] = strtok(remain,'.');
    if (isempty(remain))
        step = sscanf(token,'%d');
        break;
    end;
    basename = strcat(basename,token,'.');
end;
