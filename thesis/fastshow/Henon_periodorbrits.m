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
for n=2:12
    X=[];
    [A,num]=symbolperiod(n);
    P{1}=[P{1};num];
    for i=1:num
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
            if norm(delta,2)<1e-15 || sum>1000
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
        end
    end
    P{n}=X;
end
%save(['./data/Henon_period_orbrits_P_' num2str(x1(1)) '_' num2str(x1(2)) '_' num2str(x2(1)) '_' num2str(x2(2)) '.mat'],'P');
save(['./data/Henon_period_orbrits.mat'],'P');
%load(['./data/Henon_period_orbrits_P_' num2str(x1(1)) '_' num2str(x1(2)) '_' num2str(x2(1)) '_' num2str(x2(2)) '.mat']);
figure
C=[];
for i=4:12
    subplot(3,3,i-3)
    %figure
    count=size(P{1,i},1);
    C=[C;count];
    title(['T=' num2str(i) ', symbol=' num2str(P{1}(i)) ', num=' num2str(count)])
    hold on
    s=hsv(count);
    for j=1:count
        plot(P{1,i}(j,mod(1:end,2)==1),P{1,i}(j,mod(1:end,2)==0),'*','color',s(j,:))
    end
end
suptitle('Periodic Orbits of Henon Map')
set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
% saveas(gcf,'./temp/Henon_periodic_orbits.png')

function [A,num]=symbolperiod(n)
% 给出一个n,输出不同的符号周期轨道及符号周期轨道的个数
A=[];
for i=1:2^n-2
    A=[A;bitget(i,n:-1:1)];
end

row=1;
while row<size(A,1)
    a=A(row,:);
    for j=1:n
        a=circshift(a',1)';
        delete=[];
        for k=row+1:size(A,1)
            if a==A(k,:)
                delete=[delete,k];
            end
        end
        for l=length(delete):-1:1
            A(delete,:)=[];
        end
    end
    row=row+1;
end

b=1:n;
c=b(mod(n,b)==0);
c([1,end])=[];

row=1;
while row<=size(A,1)
    for i=1:length(c)
        times=n/c(i);
        flag=0;
        for j=2:times
            if sum(A(row,1:c(i))~=A(row,(j-1)*c(i)+1:j*c(i)))>0
                flag=1;
            end
        end
        if flag==0
            A(row,:)=[];
            break;
        end
    end
    row=row+1;
end
num=size(A,1);
end