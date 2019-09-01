function X=D1_x(n,a,q)
f=@(x)a*x;
x0=linspace(-2,2,n);
x=x0(:)';X=x0(:)';

for i=1:q
     x_temp=f(x);
     x=x_temp;
     X=[X;x];    
end
end