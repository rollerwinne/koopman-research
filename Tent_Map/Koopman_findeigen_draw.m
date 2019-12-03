function [X_eigen,X_para]=Koopman_findeigen_draw(K,U,X_eigen,X_para,n,lambda,miu,nu)
% f=@(x)norm(A*x(1:end-1)-x(end)*x(1:end-1),2)+miu*norm(diff(abs(x(1:end-1))),2);
global x
f=@(x)norm(U*x-lambda*x,2)+miu*abs(corr2(X_eigen,K*x))+nu*norm(abs(diff(x,2)));
nonlcon=@con;
options=optimoptions(@fmincon,'MaxFunEvals',10000,'Algorithm','interior-point');%'sqp');%'sqp-legacy');%
x00=zeros(length(X_para)*2,1);
x00(1:2:end)=X_para;
x00(2:2:end)=X_para;
[x,fval,exitflag,~] = fmincon(f,x00,[],[],[],[],[],[],nonlcon,options);
x_base=linspace(0,1,n);
%plot(x_base,X)
% hold on
plot(x_base,K*x)
% legend('ะก','ด๓')
title(['\lambda=',num2str(lambda),';\mu=',num2str(miu),';\nu=',num2str(nu),';fval=',num2str(fval),';exitflag=',num2str(exitflag)]);
end

function [c,ceq] = con(x)
c=[];
%c(1) = x(end)+0.1;
%c(2) = -x(end);
ceq(1) = norm(x,2)-1;
%ceq(2) = x(end);
end