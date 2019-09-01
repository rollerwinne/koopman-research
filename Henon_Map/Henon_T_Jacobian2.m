clear all;clc;%close all
q=1;a=1.4;b=0.3;
load('.\data\Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入
T=7;dot=2;
choose=1;
x0=P{T}(choose,mod(1:end,2)==1);
y0=P{T}(choose,mod(1:end,2)==0);
t=1e-8;t_temp=1e-8;
for i=1:9
    t_temp=t_temp*5;
    t=[t,t_temp]; 
end

f=@(x,y)y+1-a*x.*x;
g=@(x,y)b*x;
f_inv=@(x,y)y/b;
g_inv=@(x,y)x-1+a/b/b*y.*y;
dfgxy=@(x,y)[-2*a*x,1;b,0];

s=jet(length(x0));
str=['.\temp\test.fig'];
uiopen(str,1);  % 本征函数图像
colormap(gray)
hold on
plot3(x0,y0,0.021*ones(length(x0),1),'r*')

load(['.\data\Henon_victor_field_T4.mat']);
hold on
for i=1:length(x1)
    quiver3(x1(i),y1(i),0.021,DX1(i),DY1(i),0,'color',color1{i},'AutoScaleFactor',0.2);
    quiver3(x1(i),y1(i),0.021,-DX1(i),-DY1(i),0,'color',color1{i},'AutoScaleFactor',0.2);
    quiver3(x1(i),y1(i),0.021,DX2(i),DY2(i),0,'color',color2{i},'AutoScaleFactor',0.2);
    quiver3(x1(i),y1(i),0.021,-DX2(i),-DY2(i),0,'color',color2{i},'AutoScaleFactor',0.2);
end
axis equal

for i=1:length(x0)
    D=dfgxy(x0(i),y0(i));
    for j=i+1:i+T-1
        jj=mod(j-1,T)+1;
        D=dfgxy(x0(jj),y0(jj))*D;
    end
    [F,L]=eig(D);
    F=real(F);
    L=diag(L);
    if abs(L(1))<1
        v=F(:,1);
    elseif abs(L(2))<1
        v=F(:,2);
    else
        v=[0;0];
        disp('Error:v is null');
    end
    for k=[1,-1] %向两个方向分别演化
        for m=1:length(t) %移动不同的距离
            for n=1:dot %每个距离演化不同次数
                [X,Y]=Henon_fg_n(x0(i)+k*t(m)*v(1),y0(i)+k*t(m)*v(2),T*n,1);
                hold on
                plot3(X,Y,0.02*ones(length(X),1),'*','color',s(i,:))
            end
        end
%         [X,Y]=Henon_fg_n(x0(i)+k*t*v(1),y0(i)+k*t*v(2),T,1)
%         [X,Y]=Henon_fg_iteration(x0(i),y0(i),T,k*v(1),k*v(2),dot,0.8);

%         xx=x0(i)+k*t*v(1);
%         yy=y0(i)+k*t*v(2);
%         X=[xx];
%         Y=[yy];
%         for j=1:dot*T %演化的次数
%             x_temp=f_inv(xx,yy);
%             y_temp=g_inv(xx,yy);
%             xx=x_temp;
%             yy=y_temp;
%             X=[X,xx];
%             Y=[Y,yy];
%         end
%         hold on
%         plot3(X(1:T:end),Y(1:T:end),0.02*ones(length(X(1:T:end)),1),'linewidth',2,'color',s(i,:))
%         plot3(X,Y,0.02*ones(length(X),1),'linewidth',2,'color','green')
    end
    axis([-1.5 1.5 -1.5 1.5])
end