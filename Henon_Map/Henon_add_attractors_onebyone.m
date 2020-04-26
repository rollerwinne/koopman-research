clear;clc;close all
load('./data/Henon_attractors_data_xy.mat'); % 吸引子数据载入
load('./data/Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入
i_fig=[6,9,10,11,17,18,19];
t=[12,10,7,12,7,4,7];
choose=1;
% i_fig=6;
% T=12;

for k=1:length(t)
    str=['./temp/Henon_eigenfunctions_complex_figure' num2str(i_fig(k)) '.fig'];
    uiopen(str,1);  % 本征函数图像
    for i=1:4
        subplot(2,2,i)
        hold on
        plot3(x,y,0.02*ones(1,length(x)),'.','color',[1 0 0.8],'Markersize',1)  %吸引子图像
        hh=plot3(P{1,t(k)}(choose,mod(1:end,2)==1),P{1,t(k)}(choose,mod(1:end,2)==0),0.021*ones(1,length(P{1,t(k)}(1,:))/2),'o','color','black','Markersize',4,'MarkerFaceColor','black');
        %count=size(P{1,T},1);
        %s=jet(2);
%         for j=1:2  %周期点
%             hold on
%             if j==1
%                 hh=plot3(P{1,t(k)}(j,mod(1:end,2)==1),P{1,t(k)}(j,mod(1:end,2)==0),0.021*ones(1,length(P{1,t(k)}(1,:))/2),'o','color','black','Markersize',4,'MarkerFaceColor','black');
%             elseif j==2
%                 if size(P{1,t(k)},1)<1
%                 else
%                     % hh=plot3(P{1,T(k)}(j,mod(1:end,2)==1),P{1,T(k)}(j,mod(1:end,2)==0),0.021*ones(1,length(P{1,T(k)}(1,:))/2),'ko','Markersize',4,'MarkerFaceColor','black');
%                     % hh=plot3(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),ones(1,length(P{1,choose}(1,:))/2),'o','color','black','Markersize',4,'MarkerFaceColor','black');
%                     % hh=scatter3(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),ones(1,length(P{1,choose}(1,:))/2));
%                 end
%             end
%         end
        load(['./data/Henon_victor_field_T' num2str(t(k)) '.mat']);
        hold on
        for i=1:length(x1)   %周期点对应的本征向量
%             if strcmp(color1{i},'blue'); color1{i}='white'; end
%             if strcmp(color1{i},'red'); color1{i}='black'; end
%             if strcmp(color2{i},'blue'); color2{i}='white'; end
%             if strcmp(color2{i},'red'); color2{i}='black'; end
            quiver3(x1(i),y1(i),0.0205,DX1(i),DY1(i),0,'color',color1{i},'AutoScaleFactor',0.3);
            quiver3(x1(i),y1(i),0.0205,-DX1(i),-DY1(i),0,'color',color1{i},'AutoScaleFactor',0.3);
            quiver3(x1(i),y1(i),0.0205,DX2(i),DY2(i),0,'color',color2{i},'AutoScaleFactor',0.3);
            quiver3(x1(i),y1(i),0.0205,-DX2(i),-DY2(i),0,'color',color2{i},'AutoScaleFactor',0.3);
        end
        colormap hsv
        axis equal
    end
    drawnow
    str=['./temp/Henon_eigenfunctions_attractors_complex_figure' num2str(i_fig(k)) '.png'];
    saveas(hh,str)
end
