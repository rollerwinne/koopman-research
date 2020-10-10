clc;clear;close all;
load('Henon3d_cms_attractor.mat');
myfigure;
hold on
data=[];
index=1;
for i=1:length(Symbols)
    data1=Symbols{1,i};
    data2=Symbols{2,i};
    if isempty(data1)
        data=[data2];
    elseif isempty(data2)
        data=[data1];
    else
        data=[data1,data2];
    end
    plot3(data(1,:),data(2,:),data(3,:),'b');
    drawnow;
    pause(0.1)
end