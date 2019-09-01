close all
I=1;
deep=5;
A=abs(F(:,I));
figure
set(gcf,'outerposition',get(0,'screensize'));
plot(x0,A)
[pks,locs] = findpeaks(-A);
%locs=locs/n;
hold on
%plot(locs,A(locs),'*')
binarray=[1,2,4,8,16,32,64,128,256];

A_peak=A(locs);
[pks_temp,index]=sort(A_peak);
locs_temp=locs(index);
locs_temp=locs_temp/n;
s=jet(deep);
for i=1:deep
    which=binarray(i):binarray(i+1)-1;
    B{i}=pks_temp(which);
    B_index{i}=locs_temp(which);
    hh=plot(B_index{i},B{i},'--s',...%'LineWidth',2,...
        'MarkerSize',10,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor',s(i,:));
    for j=1:length(B{i})
        text(B_index{i}(j),B{i}(j),['(',num2str(B_index{i}(j),'%5.4f'),',',num2str(B{i}(j),'%5.4f'),')'],'VerticalAlignment','bottom','rotation',45);
    end
end

load('data0.mat');
hold on
choose=6;y_temp=0.012;
plot(X{choose},y_temp,'ko','MarkerSize',8,'MarkerFaceColor','k');
for i=1:length(X{choose})
    text(X{choose}(i),y_temp,[num2str(X{choose}(i),'%5.4f')],'VerticalAlignment','bottom','rotation',45);
end
d_abs=abs(D(h(I)));
d_angle=angle(D(h(I)))/pi*180;
str1=['n=',num2str(n),'; m=',num2str(m)];
str2=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
title({str1;str2});
str=['.\temp\Logistic_findpeaks','_n',num2str(n),'m',num2str(m),'_figure',num2str(I)];
saveas(hh,[str,'.fig'])
saveas(hh,[str,'.png'])