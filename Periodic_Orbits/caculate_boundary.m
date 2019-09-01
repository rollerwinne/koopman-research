%clear
p=6;
warning off
alpha=4;
X{1}=P{2};
for i=2:p
    X{i}=[];
    for j=1:length(X{i-1})
%         str=['4*x*(1-x)=' num2str(double(X{i-1}(1,j)))];
        str=['1-2*abs(x-1/2)=' num2str(double(X{i-1}(1,j)))];
        y=double(solve(str,'x'));
        X{i}=[X{i},y'];
    end
end
save('tent_PP.mat','X')
%save('data0.mat','X')
%save('PP.mat','X')

hold on
s=jet(p);
choose=5;y_temp=0.012;
plot(X{choose},y_temp,'ko','MarkerSize',8,'MarkerFaceColor','k');
for i=1:length(X{choose})
    text(X{choose}(i),y_temp,[num2str(X{choose}(i),'%5.4f')],'VerticalAlignment','bottom','rotation',45);
end
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
