function [R,S,T]=Henon_boundary_draw(iter,zmin,zmax,color,ab,txt)
if ab
    a=1.0;  
    b=0.54;
    B=load('Henon_boundary_ab2.mat');
else
    a=1.4;  
    b=0.3;
    B=load('Henon_boundary.mat');
end

A=B.boundary';


f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
f_inv=@(x,y)y/b;
g_inv=@(x,y)x-1+a/b/b*y.*y;

x=A(:,1);y=A(:,2);
S=[x,y];T=S;
if iter>0
    for j=1:iter
        x_temp=f(x,y);
        y_temp=g(x,y);
        x=x_temp;
        y=y_temp;
        T=[T;x,y];
    end
elseif iter<0
    for j=1:-iter
        x_temp=f_inv(x,y);
        y_temp=g_inv(x,y);
        x=x_temp;
        y=y_temp;
        T=[T;x,y];
    end
end
R=[x,y];

hold on
plot3(R(:,1),R(:,2),zmin*ones(length(R(:,1)),1),'o','Color',color,'MarkerFaceColor',color);
for i=1:length(A(:,1))
    plot3([R(i,1),R(i,1)],[R(i,2),R(i,2)],[zmin,zmax],'Color',color,'LineWidth',1);
    if nargin>5
        text(R(i,1),R(i,2),zmin,[num2str(-iter),' '],'HorizontalAlignment','right','FontSize',20);
    end
end
end