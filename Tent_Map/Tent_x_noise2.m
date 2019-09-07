%function [rou,X]=Tent_x_noise2(n,p,D)
% 测试数据
clear;close all
n=1000;%粒子个数
m=100;%基函数个数 6 100;100 10
times=10;%一次演化波包数

x0=linspace(0,1,n);
D=0.001;
f=@(x)awgn(1-2*abs(x-1/2),10*log10(1/D));%Tent map (noise)
% f=@(x)1-2*abs(x-1/2); %Tent map(clean)
% f=@(x)awgn(4.*x.*(1-x),10*log10(1/D));%Logistic map(noise)
% f=@(x)4.*x.*(1-x);%Logistic map(clean)
% figure;plot(x0,f(x0))
X=zeros(n,times);

K=zeros(n,m);
L=zeros(n,m);
for j=1:m %对于每个基函数
%     g_temp=Rect_fun(j,m);%rect基函数
    g_temp=Gauss_fun(m);%Gauss基函数
    for i=1:length(x0) %对于每个相空间的点
        K(i,j)=g_temp(x0(i),j,m);
        x_temp=[];
        for k=1:times %每个点演化迭代times次
            x_temp=[x_temp,g_temp(f(x0(i)),j,m)];
        end
        %x_temp=x_temp(x_temp<=1 & x_temp>=0);
        L(i,j)=sum(x_temp)/times;
    end
end
U=pinv(K)*L;
U=U;
[F,D]=eig(U);
D=diag(D);
% h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6);
h=1:length(D)
for i=1:min(length(h),9)
    A=abs(K*F(:,h(i)));
    figure(1);
    subplot(3,3,i)
    %set(gcf,'outerposition',get(0,'screensize'));
    hh=plot(x0,A);
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    str1=['n=',num2str(n),'; m=',num2str(m)];
    str2=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
    title({str1;str2});
end
figure(2);
SVG_draw(U,K,n,m,0,0.1,1)
figure(3);
subplot(131);spyl(U);colorbar;title('U')
subplot(132);spyl(K);colorbar;title('K')
subplot(133);spyl(L);colorbar;title('L')%K与L的矩阵形式

%
% for i=1:length(x0)
%     g_temp=Tent_Rect_fun(i,m);
%     for j=1:times
%         X(i,j)=f(x0(i));
%     end
% end
% X=X(:);
% X=X(X<=1 & X>=0);
% rou=zeros(1,n);
% roux=0:1/n:1-1/n;
% for i=1:n
%     rou(i)=sum((X>=roux(i))&(X<roux(i)+1/n));
% end
% rou(n)=rou(n)+sum(X==1);
% % rou=rou/times;
% end

function g = Rect_fun(i,m)
% x is a number
if(i==m)
    g=@(x,i,m)(x>=(i-1)/m & x<=i/m)*sqrt(m); %x坐标，第i个基函数
else
    g=@(x,i,m)(x>=(i-1)/m & x<i/m)*sqrt(m);
end
end

function g= Gauss_fun(m)
dj=1/m/2;
xj=linspace(1/2/m,1-1/2/m,m);
g=@(x,i,m)exp(-(x-xj(i)).^2./(2.*dj.^2));
end