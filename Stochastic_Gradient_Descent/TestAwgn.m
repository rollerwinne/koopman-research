clear;
n=1000;p=1;D=0.001;
x1=Tent_x(n,p);%Tent”≥…‰
x2=Tent_x_noise(n,p,D);%Tent”≥…‰ ”–‘Î…˘
%x1=Logistic_x(4,n,p);%Logistic”≥…‰
%x2=Logistic_x_noise(4,n,p,D);%Logistic”≥…‰ ”–‘Î…˘
subplot(121)
plot(x1(1,:),x1(2,:))
title('D=0')
subplot(122)
plot(x2(1,:),x2(2,:))
axis([0 1 0 1])
title(['D=',num2str(D)])

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

function X=Logistic_x_noise(alpha,n,p,D)
x0=linspace(0,1,n);
f=@(x)awgn(alpha.*x.*(1-x),10*log10(1/D));
X=x0;
x=x0;
for i=1:p
    x=f(x);
    X=[X;x];
end
end