% clear;clc;close all;
% data=load('./temp/data.mat');
% x0=data.x0;
% A=data.A;
% res=Henon_findpeaks_draw2(x0,A,'red',10);
% 
function R=Henon_findpeaks_draw(X,Y,Z,color,nearby)
%scatter3(X,Y,real(Z),3,real(Z));
%F = scatteredInterpolant(X,Y,Z);
distance=@(x,y)sqrt(sum((x-y).^2,2));
res=[];
for i=1:length(X)
    d=distance([X(i),Y(i),Z(i)],[X,Y,Z]);
    [~,idx]=sort(d);%[1,924,586,118,783,740,917,910,998,486]
    [~,index1]=sort(Z(idx(1:nearby)));%[10,8,9,7,6,2,4,3,5,1]
    [~,index2]=sort(Z(idx(1:nearby)),'descend');
    s1=idx(index1);
    s2=idx(index2);
    if s1(1)==i || s2(1)==i
        res=[res;i];
    end
end
hold on;
plot3(X(res),Y(res),Z(res),'o','Color',color,'MarkerFaceColor',color);
R=[X(res),Y(res),Z(res)];
end

