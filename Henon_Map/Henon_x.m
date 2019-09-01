function [X,Y]=Henon_x(n,a,b,q)
%% project： Henon map
tic
[x0,y0]=meshgrid(linspace(-1.5,1.5,n));
x=x0(:)';  
y=y0(:)';
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
% f=@(x,y)1-a.*abs(x)+b*y;
% g=@(x,y)x;
% f=@(x,y)[x+y+k*sin(2*pi*x),y+k*sin(2*pi*x)];
X=[x];
Y=[y];
%% 求取q次迭代
for i=1:q
     xx=f(x,y);
     yy=g(x,y); 
     x=xx;
     y=yy;
     X=[X;x]; 
     Y=[Y;y];      
end
t=toc;
disp(['Data has already been caculated and cost ' num2str(t) ' seconds.'])
%%  画演化流形图
% for i=size(X,1):-1:1
%      figure
%      plot(X(i,:),Y(i,:),'.b')
%      xlim([0,1]);
%      ylim([0,1]);
%      title(i);
% end

%save('XY.mat','X','Y');

end