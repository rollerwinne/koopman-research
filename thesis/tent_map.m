clear
alpha=2;n=1001;iter=100;
%f=@(x)x/alpha.*(x<alpha)+(1-x)./(1-alpha).*(x>alpha);
%f=@(x)x/alpha.*(x<alpha)+alpha.*(1-x).*(x>alpha);
f=@(x)alpha*min(x,1-x);
x0=linspace(0,1,n);
plot(x0,f(x0))
% n=1000;iter=100;
% for alpha=0.01:0.1:1-0.01
%     f=@(x)x/alpha.*(x<alpha)+(1-x)./(1-alpha).*(x>alpha);
%     x0=linspace(0,1,n);
%     for i=1:iter
%         x0=f(x0);
%     end
%     hold on
%     plot(alpha,x0,'r.');
% end