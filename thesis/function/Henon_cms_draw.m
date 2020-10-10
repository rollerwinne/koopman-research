function Henon_cms_draw(color,f,g,h)
X=load('Henon_cms.mat');
X=X.C0_m;
hold on
plot3(X(1,:),X(2,:),X(3,:),'o','Color',color,'MarkerFaceColor',color);
if nargin>1
    x=f(X(1,:),X(2,:),X(3,:));
    y=g(X(1,:),X(2,:),X(3,:));
    z=h(X(1,:),X(2,:),X(3,:));
    plot3(x,y,z,'o','Color',color,'MarkerFaceColor',color);
end
end