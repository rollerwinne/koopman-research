clear;clc;close all
uiopen('.\fig\Henon_eigenfunctions_complex_n100m50md45a1.4b0.3.fig',1); % 显示本征函数图像
load('.\data\Henon_attractors_data_xy.mat'); % 吸引子数据载入
load('.\data\Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入

% 绘制吸引子
for i=1:9 
    subplot(3,3,i)
    hold on
    plot(x,y,'.','color',[1 0 0.8],'Markersize',1) 
end

% 绘制周期轨道
for choose=6
    for i=1:9
        subplot(3,3,i)
%         figure
        count=size(P{1,choose},1);
        s=jet(count);
        for j=1:count
            hold on
%             hh=plot(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),'.','color','black','Markersize',15);
            if j==count
                hh=plot(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),'o','color','white','Markersize',4,'MarkerFaceColor','white');
            else
                hh=plot(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),'o','color','black','Markersize',4,'MarkerFaceColor','black');
            end
            %hh=plot(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),'o','color','black','Markersize',5);
            %axis([-1.5 1.5 -1.5 1.5])
        end
    end
end
% 
% str1=['.\fig\Henon_eigenfunctions_attractors_periodorbrits.fig'];
% str2=['.\png\Henon_eigenfunctions_attractors_periodorbrits.png'];
% saveas(hh,str1);
% saveas(hh,str2);