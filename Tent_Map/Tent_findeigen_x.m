clear;close all
n=1000;p=1;m=4;
xx=Tent_x(n,p);
% xx=Logistic_x(4,n,p);
[~,~,U,K,~,x_marker] = Tent_U_Rectangle_leftU(xx,m);
% [~,~,U,K,~,x_marker] = Tent_U_Gauss_leftU(xx,m);

A=U';
% A=[.5,.5,0,0;0,0,.5,.5;0,0,.5,.5;.5,.5,0,0];
[mm,nn]=size(A);
lambda=0;
miu=0.1;N=1;
f=@(x)norm(A*x-lambda*x)+miu*norm(diff(x),N);
nonlcon=@con;
options=optimoptions(@fmincon,'MaxFunEvals',10000)%,'Algorithm','sqp');
[x,fval,exitflag,output] = fmincon(f,2*rand(mm,1),[],[],[],[],[],[],nonlcon,options);
N1=norm(A*x-lambda*x,2)
% N2=norm(A*[.5;.5;-.5;-.5]-lambda*x,2)
D1=norm(diff(x),N)
% D2=norm(diff([.5,.5,-.5,-.5]),pp)
F1=f(x)
% F2=f([.5;.5;-.5;-.5])

pA=real(K*x);
figure
%subplot(121)
%set(gcf,'outerposition',get(0,'screensize'));
x0=linspace(0,1,n);
hh=plot(x0,pA);
hold on
plot([0,1],[0,0],'r')
title(['\mu=',num2str(miu),',\lambda=0,fval=',num2str(fval),',flag=',num2str(exitflag)])
% 
% f2=@(xx)norm(A*A*xx-lambda*xx)+miu*norm(abs(diff(xx)),N);
% nonlcon2=@con;
% [x2,fval2,exitflag2,output2] = fmincon(f2,2*rand(mm,1),[],[],[],[],[],[],nonlcon2,options);
% 
% pA2=real(K*x2);
% subplot(122)
% %set(gcf,'outerposition',get(0,'screensize'));
% x0=linspace(0,1,n);
% hh=plot(x0,pA2);
% hold on
% plot([0,1],[0,0],'r')
% title(['\mu=',num2str(miu),',\lambda=0,fval=',num2str(fval2),',flag=',num2str(exitflag2)])



function [c,ceq] = con(x)
    c=[];
    %c(1) = x(end)+0.1;
    %c(2) = -x(end);
    ceq(1) = norm(x,2)-1;
    %ceq(2) = x(end);
end