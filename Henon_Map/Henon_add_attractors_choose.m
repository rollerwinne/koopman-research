clear;clc;%close all
load('./data/Henon_attractors_data_xy.mat'); % 吸引子数据载入
load('./data/Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入
T=4;
for j=1:size(P{1,T},2)
    str=['./temp/test.fig'];
    uiopen(str,1);  % 本征函数图像
    %for i=1:4
        subplot(2,2,i)
        hold on
        plot3(x,y,0.02*ones(1,length(x)),'.','color',[1 0 0.8],'Markersize',1)  %吸引子图像
        hh=plot3(P{1,T}(j,mod(1:end,2)==1),P{1,T}(j,mod(1:end,2)==0),0.021*ones(1,length(P{1,T}(1,:))/2),'o','color','black','Markersize',4,'MarkerFaceColor','black');
        load(['./data/Henon_victor_field_T' num2str(T) '.mat']);
        hold on
        for i=1:length(x1)   %周期点对应的本征向量
            quiver3(x1(i),y1(i),0.0205,DX1(i),DY1(i),0,'color',color1{i},'AutoScaleFactor',0.3);
            quiver3(x1(i),y1(i),0.0205,-DX1(i),-DY1(i),0,'color',color1{i},'AutoScaleFactor',0.3);
            quiver3(x1(i),y1(i),0.0205,DX2(i),DY2(i),0,'color',color2{i},'AutoScaleFactor',0.3);
            quiver3(x1(i),y1(i),0.0205,-DX2(i),-DY2(i),0,'color',color2{i},'AutoScaleFactor',0.3);
        end
        colormap hsv
        axis equal
    %end
    drawnow
end