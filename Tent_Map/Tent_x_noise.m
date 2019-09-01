function X=Tent_x_noise(n,p,D)
x0=linspace(0,1,n);  
f=@(x)awgn(1-2*abs(x-1/2),10*log10(1/D));
% f=@(x)abs(1-3*abs(x-1/3)); 
X=x0;
x=x0;
for i=1:p
    x=f(x);
    X=[X;x];
end
end