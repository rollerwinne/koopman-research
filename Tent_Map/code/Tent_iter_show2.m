clear; close all;
n=100;steps=60;
f=@(x)1-2*abs(x-1/2);
x0=linspace(0,1,n);
X=x0;
for i=1:steps
    x0=f(x0);
    X=[X;x0];
end
figure
hold on
for i=1:n
    plot(1:steps+1,X(:,i))
end
