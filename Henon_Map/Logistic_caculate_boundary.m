clear
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


s=jet(p);
for j=1:9
    %subplot(3,3,j)
    hold on
%     for i=1:p
        plot(X{j},j,'*','color','b')
%     r
end

for i=1:4
    subplot(2,2,i)
    hold on
    plot(X{6},0,'*','color','r')
end
