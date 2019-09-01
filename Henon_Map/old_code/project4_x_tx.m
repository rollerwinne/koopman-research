% function [X,Y]=project4_x(n,k,q)
%% project1： 2-d map
tic
% global n;       % n:演化格点数
% global q;       % q:演化步数
% global k;       % k:非线性项
n=800;
k=0.75;
q=1000;  %迭代次数
[x0,y0]=meshgrid(linspace(0,1,n));
x=x0(:)';  
y=y0(:)';

f=@(x,y)mod((x+y+k/(2*pi)*sin(2*pi*x)),1);
g=@(x,y)mod((y+k/(2*pi)*sin(2*pi*x)),1);
% f=@(x,y)[x+y+k*sin(2*pi*x),y+k*sin(2*pi*x)];
X=[x];
Y=[y];
%% 求取q次迭代
for i=1:q 
%      y=mod(g(x,y),1);
%      x=mod(f(x,y),1);

     xx=f(x,y);
     yy=g(x,y); 
     x=xx;
     y=yy;
     X=[X;x]; 
     Y=[Y;y];      
end

%%  画演化流形图
% for i=size(X,1):-1:1
%      figure
%      plot(X(i,:),Y(i,:),'.b')
%      xlim([0,1]);
%      ylim([0,1]);
%      title(i);
% end
save('XY.mat','X','Y');
%% 本征函数为调谐时间平均值函数 （本征值为1）
F=0;
for i=1:q
   h=cos(2*pi.*X(i,:)).*cos(2*pi.*Y(i,:));
    F=F+h;
end
    F=F/q;
%% 画图
    figure('NumberTitle','off','Name',['q=' num2str(q) 'k=' num2str(k) 'n=' num2str(n)]);
    for i=1:q
    hh=scatter3(X(i,:)',Y(i,:)',F',3,F');
    colorbar
    colormap(jet)
    view(0,90)
    axis equal
    end
    str=['n' num2str(n) 'q' num2str(q)  'k' num2str(k) '.fig'];
    saveas(hh,str);
    
    toc
% end
