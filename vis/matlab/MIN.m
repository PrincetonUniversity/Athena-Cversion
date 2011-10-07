function y = MIN(x)
y = x;
for i = 1:ndims(x)
    y = min(y);
end;