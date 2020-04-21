clear
iter=10000;
for alpha=0.001:0.001:1-0.001
    f=@(x)x/alpha.*(x<alpha)+(1-x)./(1-alpha).*(x>alpha);
    df=@(x)1/alpha.*(x<alpha)-1/(1-alpha).*(x>alpha);
    s=zeros(1,iter);
    x0=0.2;
    for i=1:iter
        x0=f(x0);
        s(i)=df(x0);
    end
    lypn=1/iter*sum(log(s));
    hold on
    plot(alpha,lypn,'b.')
end
plot([0,1],[0,0],'r')
ylabel('\lambda')
xlabel('\alpha')
title('Lyapunov Exponent of Tent Map')
