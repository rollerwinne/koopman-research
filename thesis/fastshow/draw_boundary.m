clear
load('logistic_boundary_norepeat_x0.mat', 'X')
for i=1:length(X)
    Y{i}=sort(X{i});
end
draw1D(Y,1,5);
title('Boundary of Logistic Map (x=0)')
saveas(gcf,'./temp/Logistic_boundarys_x0.png')

function draw1D(X,m,n)
yshift=0.05;
hold on
plot([0 1],[0 0],'k-o','MarkerSize',1,'MarkerFaceColor','k');
text(0,0-yshift,'0','Horiz','center','Vert','top','FontSize',20);
text(1,0-yshift,'1','Horiz','center','Vert','top','FontSize',20);
drawBoundary(X,m,n)
axis off
end

function drawBoundary(X,m,n)
yshift=0.05;
size=15;rate=0.87;
s=flipud(jet(n-m+1));
for j=m:n
    x=X{j};
    plot(x,0,'bo','MarkerSize',size*rate^j,'MarkerFaceColor',s(j,:))
    for i=1:length(x)
        text(x(i),0+yshift,num2str(x(i),'%.4f'),'rotation',45)
        text(x(i),0-yshift,num2str(j));
    end
end
end