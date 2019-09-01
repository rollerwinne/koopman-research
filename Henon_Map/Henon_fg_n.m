function [x,y]=Henon_fg_n(x,y,n,p)
a=1.4;b=0.3;
f=@(x,y)y+1-a*x.*x;
g=@(x,y)b*x;
f_inv=@(x,y)y/b;
g_inv=@(x,y)x-1+a/b/b*y.*y;
if p==0 %正向迭代
    for i=1:n
        x_temp=f(x,y);
        y_temp=g(x,y);
        x=x_temp;
        y=y_temp;
    end
elseif p==1 %反向迭代
    for i=1:n
        x_temp=f_inv(x,y);
        y_temp=g_inv(x,y);
        x=x_temp;
        y=y_temp;
    end
end
end