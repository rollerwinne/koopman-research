function X=Tent_x(n,p)
x0=linspace(0,1,n);  
f=@(x)1-2*abs(x-1/2);
% f=@(x)abs(1-3*abs(x-1/3)); 
X=x0;
x=x0;
for i=1:p
    x=f(x);
    X=[X;x];
end
end