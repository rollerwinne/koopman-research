clear;close all
n=1000;p=1;mi=1;
for m=[4,6,7,10,15,20,30,50,60,75,100,200]
    %xx=Tent_x(n,p);%Tent”≥…‰
    %xx=Tent_x_noise(n,p,0.001);%Tent”≥…‰ ”–‘Î…˘
    %xx=Logistic_x(4,n,p);%Logistic”≥…‰
    xx=Logistic_x_noise(4,n,p,0.001);%Logistic”≥…‰ ”–‘Î…˘
    [~,~,U,K,~,x_marker] = Tent_U_Rectangle_leftU(xx,m);%Rectangleª˘∫Ø ˝
    % [~,~,U,K,~,x_marker] = Tent_U_Gauss_leftU(xx,m);%∏ﬂÀπª˘∫Ø ˝
    
    A=U';
    % A=[.5,.5,0,0;0,0,.5,.5;0,0,.5,.5;.5,.5,0,0];
    [mm,nn]=size(A);
    lambda=0;N=1;
    miu=0.1;nu=1;
    % full_x=@(x)[x;1-norm(x,N)];
    f=@(x)norm(A*x-lambda*x)+miu*norm(diff(x),N)+nu*(abs(norm(x)-1));
    
    X=[];
    for i=1:1
        x0=rand(mm,1);yita=0.005;iter=3000;
        x=fmin_GD(f,x0,yita,iter);
        X=[X,x];
        %f(x)
    end
    x=sum(X,2)/length(X(1,:));
    
    pA=real(K*x);
    %figure
    subplot(3,4,mi);mi=mi+1;
    %set(gcf,'outerposition',get(0,'screensize'));
    x0=linspace(0,1,n);
    hh=plot(x0,pA);
    hold on
    plot([0,1],[0,0],'r')
    title(['m=',num2str(m),',fval=',num2str(f(x))])
end
suptitle(['p=',num2str(N),',\lambda=0,\mu=',num2str(miu),',\nu=',num2str(nu)])

function x=fmin_GD(f,x0,yita,iter)
x=x0;
%disp(x);
for i=1:iter
    dx=zeros(length(x0),1);
    for j=1:length(x0)
        dx(j)=DFx(f,x,j,0.001);
    end
    x=x-yita*dx;
end
end

function dx=DFx(f,x,i,delta_x)
dx=(f(x+[zeros(i-1,1);delta_x;zeros(length(x)-i,1)])-f(x))/delta_x;
end

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

function [F,D,U,K,L,x0]= Tent_U_Rectangle_leftU(x,m,flag)
x_k=x(1,:);
x_l=x(end,:);
x0=linspace(1/2/m,1-1/2/m,m);
if nargin<3
    K=Logistic_G_Rectangle(x_k,m);
    L=Logistic_G_Rectangle(x_l,m);
else
    K=Logistic_G_Rectangle_natural(x_k,m);
    L=Logistic_G_Rectangle_natural(x_l,m);
end
U=pinv(K)*L;
[F,D]=eig(U);
D=diag(D);
end

function G= Logistic_G_Rectangle(x,m)
G=zeros(length(x),m);
x=x(:);
for i=1:m-1
    G(:,i)=(x>=(i-1)/m & x<i/m)*sqrt(m);
end
G(:,m)=(x>=(m-1)/m & x<=m/m)*sqrt(m);
end

function G=Logistic_G_Rectangle_natural(x,m)
x_function=Logistic_x(4,m+1,1);
x_function(2,:);
for i=1:m-1
    G(:,i)=(x>=x_function(2,i) & x<x_function(2,i+1)) * sqrt(1/(x_function(2,i+1)-x_function(2,i)));
end
i=m;
G(:,i)=(x>=x_function(2,i) & x<x_function(2,i+1)) * sqrt(1/(x_function(2,i+1)-x_function(2,i)));
end