clear all;clc;close all;tic
q=1;a=1.4;b=0.3;
T=7;
for choose=3:4;
    load('.\data\Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入
    x0=P{T}(choose,mod(1:end,2)==1);
    y0=P{T}(choose,mod(1:end,2)==0);
    t=0.5e-7;dot=1;n=10000;
    
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
    
    load(['.\data\Henon_victor_field_T' num2str(T) '_choose' num2str(choose) '.mat']);
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
            disp('Error:cannot find stable eigenvalue');
        end
        k_d=[1,-1];
        for k=1:length(k_d) %向两个方向分别演化
            [X{i,k},Y{i,k}]=Henon_manifold_computing3(x0(i),y0(i),k_d(k)*v(1),k_d(k)*v(2),T,t,n,dot);
            time=toc;
            disp(['i=' num2str(i)  ', k=' num2str(k_d(k)) ' is done. Used ' num2str(time) ' seconds']);
            hold on
            plot3(X{i,k},Y{i,k},0.02*ones(length(X{i,k}),1),'LineWidth',2,'color',s(i,:))
            %plot3(X{i,k},Y{i,k},0.02*ones(length(X{i,k}),1),'*','color',s(i,:))
            drawnow
        end
        axis([-1.5 1.5 -1.5 1.5])
    end
    xlabel(['t=' num2str(t) '; n=' num2str(n) '; dot=' num2str(dot)])
end