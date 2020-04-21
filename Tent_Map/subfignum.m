function [x,y]=subfignum(n)
y=ceil(sqrt(n));
if y*(y-1)<n
    x=y;
else
    x=y-1;
end
end