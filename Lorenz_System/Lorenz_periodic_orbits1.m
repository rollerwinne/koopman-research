clear %;clf,close all
% unable to use
%% Lorenz equation
rho=28;sigma=10;beta=8/3;
f=@(x,y,z)[sigma*(y-x);x.*(rho-z)-y;x.*y-beta*z];
df=@(x,y,z)[-sigma,sigma,0;
    rho-z,-1,-x;
    y,x,-beta];
%% Initial point
x1=[-15,-18,34];
x2=[10,18,18];
%% Multipoint shooting
P{1}=[0];
for n=2:10
    X=[];
    [A,num]=symbolperiod(n);
    P{1}=[P{1};num];
    for i=1:num
        xx=[];
        yy=[];
        zz=[];
        for l=1:n
            if A(i,l)==1
                x=[x;x1'];
                xx=[xx;x1(1)];
                yy=[yy;x1(2)];
                zz=[zz;x1(3)];
            elseif A(i,l)==0
                x=[x;x2'];
                xx=[xx;x2(1)];
                yy=[yy;x2(2)];
                zz=[zz;x2(3)];
            end
        end
        sum=0;
        while true
            F=[];
            DF=[];
            for j=2:n
                x_temp=rk_4(f,[0,0.01,0.01],xx(j-1),yy(j-1),zz(j-1));
                F=[F;xx(j)-x_temp(end,1);yy(j)-x_temp(end,2);zz(j)-x_temp(end,3)];
                DF=[DF;[zeros(2,2*(j-2)),-dfgxy(xx(j-1),yy(j-1)),eye(2),zeros(2,2*(n-j))]];
            end
            F=[xx(1)-f(xx(n),yy(n));yy(1)-g(xx(n),yy(n));F];
            DF=[eye(2),zeros(2,2*(n-2)),-dfgxy(xx(n),yy(n));DF];
            delta=-DF\F;
            sum=sum+1;
            if norm(delta,2)<1e-15 || sum>1000;
                break;
            else
                x=x+delta;
                xx=[];
                yy=[];
                for k=1:n
                    xx=[xx;x(2*k-1)];
                    yy=[yy;x(2*k)];
                end
            end
        end
        if sum<=1000
            X=[X;x'];
            disp(['n=' num2str(n) ' num=' num2str(i) '/' num2str(num) ' ÊÕÁ²']);
        end
    end
    P{n}=X;
end
save(['.\data\Henon_period_orbrits_P_' num2str(x1(1)) '_' num2str(x1(2)) '_' num2str(x2(1)) '_' num2str(x2(2)) '.mat'],'P');
