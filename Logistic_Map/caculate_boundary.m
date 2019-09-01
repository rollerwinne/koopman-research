%clear
p=9;
warning off
alpha=4;
X{1}=0;
for i=2:p
    X{i}=[];
    for j=1:length(X{i-1})
        str=['4*x*(1-x)=' num2str(double(X{i-1}(1,j)))];
        y=double(solve(str,'x'));
        X{i}=[X{i},y'];
    end
end
save('data0.mat','X')

hold on
s=jet(p);
choose=6;y_temp=0.012;
plot(X{choose},y_temp,'ko','MarkerSize',8,'MarkerFaceColor','k');
for i=1:length(X{choose})
    text(X{choose}(i),y_temp,[num2str(X{choose}(i),'%5.4f')],'VerticalAlignment','bottom','rotation',45);
end
d_abs=abs(D(h(i)));
d_angle=angle(D(h(i)))/pi*180;
str1=['n=',num2str(n),'; m=',num2str(m)];
str2=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
title({str1;str2});
str=['.\temp\Logistic_findpeaks','_n',num2str(n),'m',num2str(m),'_figure',num2str(i)];
saveas(hh,[str,'.png'])
% for j=1:9
%     %subplot(3,3,j)
%     hold on
% %     for i=1:p
%         plot(X{j},j,'*','color','b')
% %     r
% end
% 
% for i=1:4
%     subplot(2,2,i)
%     hold on
%     plot(X{6},0,'*','color','r')
% end
