clear all;clc;close all
n=20;q=1;a=1.4;b=0.3;
%% all phase space
% [x0,y0]=meshgrid(linspace(-1.5,1.5,n));
% x1=x0(:)';
% y1=y0(:)';
%% only on period orbrits
load('.\data\Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入
for choose=1:4
    T=7;
    x1=P{T}(choose,mod(1:end,2)==1);
    y1=P{T}(choose,mod(1:end,2)==0);
    %%
    f=@(x,y)y+1-a*x.*x;
    g=@(x,y)b*x;
    % for T=9;
    syms x y
    for i=1:T
        x_temp=f(x,y);
        y_temp=g(x,y);
        x=x_temp;
        y=y_temp;
    end
    X=x;
    Y=y;
    syms x y
    J=[diff(X,x),diff(X,y);diff(Y,x),diff(Y,y)];
    
    for i=1:length(x1)
        D=double(subs(J,[x,y],[x1(i),y1(i)]));
        [F,L]=eig(D);
        F=real(F);
        L=diag(L);
        if abs(L(1))<1
            color1{i}='blue';
        else
            color1{i}='red';
        end
        if abs(L(2))<1
            color2{i}='blue';
        else
            color2{i}='red';
        end
        DX1(i)=F(1,1);
        DY1(i)=F(2,1);
        DX2(i)=F(1,2);
        DY2(i)=F(2,2);
    end
    save(['.\data\Henon_victor_field_T' num2str(T) '_choose' num2str(choose) '.mat'],'T','x1','y1','DX1','DX2','DY1','DY2','color1','color2');
    %subplot(3,3,T)
    figure(choose)
    hold on
    for i=1:length(x1)
        quiver3(x1(i),y1(i),0.021,DX1(i),DY1(i),0,'color',color1{i},'AutoScaleFactor',0.2);
        quiver3(x1(i),y1(i),0.021,-DX1(i),-DY1(i),0,'color',color1{i},'AutoScaleFactor',0.2);
        quiver3(x1(i),y1(i),0.021,DX2(i),DY2(i),0,'color',color2{i},'AutoScaleFactor',0.2);
        quiver3(x1(i),y1(i),0.021,-DX2(i),-DY2(i),0,'color',color2{i},'AutoScaleFactor',0.2);
    end
    colormap hsv
    axis equal
    %title(['p=' num2str(T)])
    %plot3(x1,y1,0.03*ones(1,length(x1)),'o','color','black','Markersize',4,'MarkerFaceColor','black');
    drawnow
    % end
end