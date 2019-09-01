function X=Logistic_x(alpha,n,p)
x0=linspace(0,1,n);  
f=@(x)alpha.*x.*(1-x);
X=x0;
x=x0;
for i=1:p
    x=f(x);
    X=[X;x];
end
end