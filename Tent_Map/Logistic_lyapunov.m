clear;
n=10000;
a=3:0.001:4;
len=length(a);
a=reshape(a,len,1);
sum=zeros(len,1);
unit=ones(len,1);
x=unit*0.1;
 
for i=1:n
    y=2-(x>0.5)*4;
    sum=sum+log(abs(y));
    x=1-2*abs(x-1/2);
end
lamuda=sum/n;
plot(a,lamuda)
grid on

xlabel('\fontsize{16}a')
ylabel('\fontsize{16}Lyapunov\fontname{隶书}指数\lambda')
title('\fontsize{16}a\fontname{隶书}指数Lyapunov\fontname{隶书}指数的关系曲线')