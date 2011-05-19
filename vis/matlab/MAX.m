function y = MAX(x)
y = x;
for i = 1:ndims(x)
    y = max(y);
end;