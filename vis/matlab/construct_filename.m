function filename = construct_filename(path,basename,step)
% 
% construct_filename:  CONSTRUCT A FULL-PATH FILENAME FROM ITS COMPONENT
% PARTS--THAT WAY, ONE NEED ONLY SPECIFY A SINGLE PATH, BASENAME, AND STEP
% RANGE TO OPEN MULTIPLE FILES.
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/22/09

filename = strcat(path,basename,'.',sprintf('%04d',step),'.bin');