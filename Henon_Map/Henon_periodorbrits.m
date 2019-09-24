clear %;clf,close all
a=1.4;
b=0.3;
%% Lozi map equation
% f=@(x,y)1-a.*abs(x)+b*y;
% g=@(x,y)x;
% dfx=@(x,y)a*(x<0)+0*(x==0)-a*(x>0);
% dfy=@(x,y)b;
% dgx=@(x,y)1;
% dgy=@(x,y)0;
% dfgxy=@(x,y)[a*(x<0)+0*(x==0)-a*(x>0),b;1,0];

%% Henon map equation
f=@(x,y)y+1-a*x.*x;
g=@(x,y)b*x;
dfx=@(x,y)-2*a*x;
dfy=@(x,y)1;
dgx=@(x,y)b;
dgy=@(x,y)0;
dfgxy=@(x,y)[-2*a*x,1;b,0];
%% Initial point
% x1=[-1,0.2];
% x2=[1,-0.2];
x1=[1,0.3];
x2=[-1,-0.3];
% x1=[1,1];
% x2=[-1,-1];
% x1=[-1,1];
% x2=[1,-1];

%%
P{1}=[0];
for n=2:10
    X=[];
    [A,num]=symbolperiod(n);
    P{1}=[P{1};num];
    for i=1:num
        %     赋随机初值（通过迭代）
        %     x=zeros(2*n,1);
        %     x(1:2)=-1.5+3*rand(2,1);
        %     xx=x(1);
        %     yy=x(2);
        %     for i=2:n
        %         x(2*i-1)=f(2*i-3,2*i-2);
        %         xx=[xx;x(2*i-1)];
        %         x(2*i)=g(2*i-3,2*i-2);
        %         yy=[yy;x(2*i)];
        %     end
        x=[];
        xx=[];
        yy=[];
        for l=1:n
            if A(i,l)==1
                x=[x;x1(1);x1(2)];
                xx=[xx;x1(1)];
                yy=[yy;x1(2)];
            elseif A(i,l)==0
                x=[x;x2(1);x2(2)];
                xx=[xx;x2(1)];
                yy=[yy;x2(2)];
            end
        end
        sum=0;
        while true
            F=[];
            DF=[];
            for j=2:n
                F=[F;xx(j)-f(xx(j-1),yy(j-1));yy(j)-g(xx(j-1),yy(j-1))];
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
            disp(['n=' num2str(n) ' num=' num2str(i) '/' num2str(num) ' 收敛']);
            %             subplot(3,3,n-1)
            %             plot(xx,yy,'*','color','r')
            %             %axis([-1.5 1.5 -1.5 1.5])
            %             title(n)
        else
            % disp(['n=' num2str(n) ' num=' num2str(i) '/' num2str(num) ' 不收敛']);
            %             subplot(3,3,n-1)
            %             plot(0,0,'x')
            %             axis([-1.5 1.5 -1.5 1.5])
            %             title(n)
        end
    end
    P{n}=X;
end
save(['./data/Henon_period_orbrits_P_' num2str(x1(1)) '_' num2str(x1(2)) '_' num2str(x2(1)) '_' num2str(x2(2)) '.mat'],'P');
%load(['./data/Henon_period_orbrits_P_' num2str(x1(1)) '_' num2str(x1(2)) '_' num2str(x2(1)) '_' num2str(x2(2)) '.mat']);
figure
C=[];
for i=2:10
    subplot(3,3,i-1)
    %figure
    count=size(P{1,i},1);
    C=[C;count];
    title(['T=' num2str(i) ',count=' num2str(count) ',Total=' num2str(P{1}(i))])
    hold on
    s=jet(count);
    for j=1:count
        %figure
        plot(P{1,i}(j,mod(1:end,2)==1),P{1,i}(j,mod(1:end,2)==0),'*','color',s(j,:))
        %axis([-1.5 1.5 -1.5 1.5])
    end
end