clear %;clf,close all
% unable to use
%% Lorenz equation
rho=28;sigma=10;beta=8/3;
f=@(x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
df=@(x)[-sigma,sigma,0;
    rho-x(3),-1,-x(1);
    x(2),x(1),-beta];
%% Initial point
x0=[-1,3,4];
x1=[12,17,27];
x2=[-13,-18,27];
% x0=[-1,3,4];
% x1=Lorenz_Poincare_next_point(x0);
%% Multipoint shooting
P{1}=[0];
for n=8 % 对每个周期循环
    X=[];
    [A,num]=symbolperiod(n);
    P{1}=[P{1};num];
    for i=1:num % 对每个周期的每个符号周期循环
        x=[];xx=[];
        for l=1:n % 设初值
            if A(i,l)==1
                x=[x;x1'];
                xx=[xx;x1];
            elseif A(i,l)==0
                x=[x;x2'];
                xx=[xx;x2];
            end
        end
        sum=0;
        while true
            F=[];
            DF=[];
            for j=2:n
                xx_temp=Lorenz_Poincare_next_point(xx(j-1,:));
                F=[F;[xx(j,:)-xx_temp]'];
                DF=[DF;[zeros(3,3*(j-2)),-df(xx(j-1,:)),eye(3),zeros(3,3*(n-j))]];
            end
            xx_temp=Lorenz_Poincare_next_point(xx(n,:));
            F=[[xx(1,:)-xx_temp]';F];
            DF=[eye(3),zeros(3,3*(n-2)),-df(xx(n,:));DF];
            delta=-DF\F;
            sum=sum+1;
            if norm(delta,2)<1e-8 || sum>100
                break;
            else
                disp(['n=',num2str(n),' num=',num2str(i),'/',num2str(num),' count=',num2str(sum)]);
                delta(mod(1:end,3)==0)=0;
                x=x+delta;
                xx(:,1)=x(mod(1:end,3)==1);
                xx(:,2)=x(mod(1:end,3)==2);
                xx(:,3)=x(mod(1:end,3)==0);
            end
        end
        if sum<=100
            X=[X;x'];
            disp(['n=' num2str(n) ' num=' num2str(i) '/' num2str(num) ' 收敛']);
        else
            disp(['n=' num2str(n) ' num=' num2str(i) '/' num2str(num) ' 不收敛']);
        end
    end
    P{n}=X;
end
%save(['.\data\Henon_period_orbrits_P_' num2str(x1(1)) '_' num2str(x1(2)) '_' num2str(x2(1)) '_' num2str(x2(2)) '.mat'],'P');
