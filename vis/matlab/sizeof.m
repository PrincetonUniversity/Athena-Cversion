function nbytes = sizeof(datatype)
% 
% sizeof:  RETURNS THE SIZE, IN BYTES, OF A GIVEN DATA TYPE.
%
% AUTHOR:  AARON SKINNER
% LAST MODIFIED:  6/29/09

try
    z = zeros(1,datatype);
catch
    error('[sizeof]:  Unsupported data type!');
end

w = whos('z');
nbytes = w.bytes;