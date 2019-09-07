function SVG_draw(U,K,n,m,lambda,miu,nu)
A=U';
[mm,nn]=size(A);
lambda=0;N=1;
miu=0.1;nu=1;
f=@(x)norm(A*x-lambda*x)+miu*norm(diff(x),N)+nu*(abs(norm(x)-1));

X=[];
for i=1:1
    x0=rand(mm,1);yita=0.005;iter=3000;
    x=fmin_GD(f,x0,yita,iter);
    X=[X,x];
end
x=sum(X,2)/length(X(1,:));

pA=real(K*x);
% figure
%subplot(3,4,mi);mi=mi+1;
%set(gcf,'outerposition',get(0,'screensize'));
x0=linspace(0,1,n);
hh=plot(x0,pA);
hold on
plot([0,1],[0,0],'r')
title(['SVG_m=',num2str(m),',fval=',num2str(f(x))])
end

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

