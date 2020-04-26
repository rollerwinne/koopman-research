function Henon_T_draw(T,choose,yt)
period=load('./data/Henon_period_orbrits.mat'); % 周期轨道数据载入
P=period.P;
count=size(P{1,T},1);
X=P{1,T}(choose,mod(1:end,2)==1);
Y=P{1,T}(choose,mod(1:end,2)==0);
Z=yt*ones(size(X));
hh=plot3(X,Y,Z,'o','color','black','Markersize',4,'MarkerFaceColor','black');
view(0,90)
end