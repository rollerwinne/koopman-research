function [rou,X]=Tent_x_noise2(n,p,D)
x0=linspace(0,1,n);
times=100;
m=100;%基函数个数
% f=@(x)awgn(1-2*abs(x-1/2),10*log10(1/D));
f=@(x)1-2*abs(x-1/2);
X=zeros(n,times);
for i=1:length(x0)
    g_temp=Tent_Rect_fun(i,m);
    for j=1:times
        X(i,j)=f(x0(i));
    end
end
X=X(:);
X=X(X<=1 & X>=0);
rou=zeros(1,n);
roux=0:1/n:1-1/n;
for i=1:n
    rou(i)=sum((X>=roux(i))&(X<roux(i)+1/n));
end
rou(n)=rou(n)+sum(X==1);
% rou=rou/times;
end

function g = Tent_Rect_fun(i,m)
% x is a number
if(i==m)
    g=@(x)(x>=(i-1)/m & x<=i/m)*sqrt(m);
else
    g=@(x)(x>=(i-1)/m & x<i/m)*sqrt(m);
end
end