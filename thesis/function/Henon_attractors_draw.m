function Henon_attractors_draw(color,zt,~)
if nargin<1
    color=[1 0 0.8];
end
attr=load('./data/Henon_attractors_data_xy.mat'); % 吸引子数据载入
hold on
plot3(attr.x,attr.y,zt*ones(size(attr.x)),'.','color',color,'Markersize',1)
if nargin>2
    plot3([attr.x1,attr.x2],[attr.y1,attr.y2],[zt,zt],'r*');
end
end