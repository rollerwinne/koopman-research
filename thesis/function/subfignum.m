function [x,y]=subfignum(n)
y=ceil(sqrt(n));
if y*(y-1)<n
    x=y;
else
    x=y-1;
end
if x>3
    x=3;
end
if y>3
    y=3;
end
end