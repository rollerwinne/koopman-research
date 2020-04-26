function [X,Y]=Henon_x(n,a,b,q)
[x0,y0]=meshgrid(linspace(-1.5,1.5,n));
x=x0(:)';  
y=y0(:)';
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
X=[x];
Y=[y];
for i=1:q
     x_temp=f(x,y);
     y_temp=g(x,y); 
     x=x_temp;
     y=y_temp;
     X=[X;x]; 
     Y=[Y;y];      
end

% figure
% for i=1:size(X,1)
%      plot(X(i,:),Y(i,:),'.b')
%      title(i);
%      axis([-1.5,1.5,-1.5,1.5]);
%      drawnow;
%      pause(0.3)
% end

end