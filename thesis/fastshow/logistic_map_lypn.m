clear
iter=10000;
for gamma=3:0.001:4
    f=@(x)gamma.*x.*(1-x);
    df=@(x)gamma-2*gamma*x;
    s=zeros(1,iter);
    x0=0.1;
    
    for i=1:iter
        x0=f(x0);
        s(i)=df(x0);
    end
    lypn=1/iter*sum(log(s));
    hold on
    plot(gamma,lypn,'b.')
end
plot([3,4],[0,0],'r')
ylabel('\lambda')
xlabel('\gamma')
title('Lyapunov Exponent of Logistic Map')
