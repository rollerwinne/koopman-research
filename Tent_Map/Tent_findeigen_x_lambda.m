clear
n=1000;p=1;m=4;
xx=Tent_x(n,p);
[~,~,U,K,~,x_marker] = Tent_U_Rectangle_leftU(xx,m);
A=U;
[mm,nn]=size(A);
miu=0.25;
f=@(x)norm(A*x(1:end-1)-x(end)*x(1:end-1),2)+miu*norm(diff(abs(x(1:end-1))),2);
nonlcon=@con;
[x,fval,exitflag,output] = fmincon(f,rand(mm+1,1),[],[],[],[],[],[],nonlcon)
% N1 = norm(A*x(1:end-1)-x(end)*x(1:end-1),2)
% D1 = norm(diff(x(1:end-1)),2)
% D2 = norm([.5,.5,-.5,-.5],2)
% f(x)

pA=real(K*x(1:end-1));
figure
%set(gcf,'outerposition',get(0,'screensize'));
x0=linspace(0,1,n);
hh=plot(x0,pA);
hold on
plot([0,1],[0,0],'r')
title(['\mu=',num2str(miu),',\lambda=',num2str(x(end)),',fval=',num2str(fval),',flag=',num2str(exitflag)])

f2=@(xx)norm(A*A*xx(1:end-1)-xx(end)*xx(1:end-1)-x(1:end-1))+miu*norm(abs(diff(xx(1:end-1),2)));
nonlcon2=@con;
[x2,fval2,exitflag2,output2] = fmincon(f2,rand(mm+1,1),[],[],[],[],[],[],nonlcon2);

pA2=real(K*x2(1:end-1));
figure
%set(gcf,'outerposition',get(0,'screensize'));
x0=linspace(0,1,n);
hh=plot(x0,pA2);
hold on
plot([0,1],[0,0],'r')
title(['\mu=',num2str(miu),',\lambda=0,fval=',num2str(fval2),',flag=',num2str(exitflag2)])

