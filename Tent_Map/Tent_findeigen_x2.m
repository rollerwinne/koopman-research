clear;close all
n=1000;p=1;mi=1;
for m=[4,6,7,10,15,20,30,40,50,60,80,100]
    %xx=Logistic_x(4,n,p);
    xx=Tent_x(n,p);
    [~,~,U,K,~,x_marker] = Tent_U_Rectangle_leftU(xx,m);
    %[~,~,U,K,~,x_marker] = Tent_U_Gauss_leftU(xx,m);
    A=U';
    
    [mm,nn]=size(A);
    lambda=0;
    miu=0.1;N=1;
    f=@(x)norm(A*x-lambda*x)+miu*norm(diff(x),N);
    nonlcon=@con;
    options=optimoptions(@fmincon,'MaxFunEvals',10000,'Algorithm','sqp-legacy');%'sqp');
    [x,fval,exitflag,output] = fmincon(f,2*rand(mm,1),[],[],[],[],[],[],nonlcon,options);
    N1=norm(A*x-lambda*x,2)
    D1=norm(diff(x),N)
    F1=f(x)
    
    subplot(3,4,mi);mi=mi+1;
    pA=real(K*x);
    x0=linspace(0,1,n);
    hh=plot(x0,pA);
    hold on
    plot([0,1],[0,0],'r')
    title(['m=',num2str(m),',fval=',num2str(fval),',flag=',num2str(exitflag)])
end
suptitle(['\mu=',num2str(miu),',\lambda=0'])

function [c,ceq] = con(x)
c=[];
%c(1) = x(end)+0.1;
%c(2) = -x(end);
ceq(1) = norm(x,2)-1;
%ceq(2) = x(end);
end