clear;clc;close all
load('.\data\Henon_attractors_data_xy.mat'); % 吸引子数据载入
load('.\data\Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入
choose=6;

for i=1:9
    str=['.\temp\Henon_eigenfunctions_complex_figure' num2str(i) '.fig'];
    uiopen(str,1);
    hold on
    plot3(x,y,0.02*ones(1,length(x)),'.','color',[1 0 0.8],'Markersize',1)
    count=size(P{1,choose},1);
    s=jet(count);
    for j=1:count
        hold on
        if j==count
        hh=plot3(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),0.021*ones(1,length(P{1,choose}(1,:))/2),'o','color','white','Markersize',4,'MarkerFaceColor','white');
        else
        hh=plot3(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),0.021*ones(1,length(P{1,choose}(1,:))/2),'ko','Markersize',4,'MarkerFaceColor','black');
        % hh=plot3(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),ones(1,length(P{1,choose}(1,:))/2),'o','color','black','Markersize',4,'MarkerFaceColor','black');
        % hh=scatter3(P{1,choose}(j,mod(1:end,2)==1),P{1,choose}(j,mod(1:end,2)==0),ones(1,length(P{1,choose}(1,:))/2));
        end
    end
    drawnow
    str=['.\temp\real_figure' num2str(i) '.png'];
    saveas(hh,str)
end

% str=['.\temp\real_figure' num2str(i) '.png'];
% saveas(hh,str)

% str1=['.\fig\Henon_eigenfunctions_attractors_periodorbrits.fig'];
% str2=['.\png\Henon_eigenfunctions_attractors_periodorbrits.png'];
% saveas(hh,str1);
% saveas(hh,str2);