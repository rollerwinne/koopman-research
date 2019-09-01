function [X,Y]=D2_x(n,a,b,q)
f=@(x,y)a*x;
g=@(x,y)b*y;
[x0,y0]=meshgrid(linspace(-2,2,n));
x=x0(:)';X=x0(:)';
y=y0(:)';Y=y0(:)';

for i=1:q
     xx=f(x,y);
     yy=g(x,y); 
     x=xx;
     y=yy;
     X=[X;x]; 
     Y=[Y;y];      
end
end