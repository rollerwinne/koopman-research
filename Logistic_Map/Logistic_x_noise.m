function X=Logistic_x_noise(alpha,n,p,e)
x0=linspace(0,1,n);  
f=@(x)alpha.*x.*(1-x);
X=x0;
x=x0;
for i=1:p
    x=f(x)+e*randn(1,n);
    X=[X;x];
end
end