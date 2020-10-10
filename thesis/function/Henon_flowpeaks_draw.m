function R=Henon_flowpeaks_draw(X,Y,Z,color)
[~,loc1]=findpeaks(Z);
[~,loc2]=findpeaks(-Z);
loc=[loc1;loc2];
hold on;
plot3(X(loc),Y(loc),Z(loc),'o','Color',color,'MarkerFaceColor',color);
plot3(X(loc),Y(loc),min(Z(:))*ones(size(X(loc))),'o','Color',color,'MarkerFaceColor',color);
R=[X(loc),Y(loc),Z(loc)];
end