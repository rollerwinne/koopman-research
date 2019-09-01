clear;
D=100;
f=@(x)sin(x);
g=@(x)awgn(sin(x),10*log10(D));
x=linspace(0,2*pi,1000);
plot(x,f(x),'b.')
hold on
plot(x,g(x),'r')