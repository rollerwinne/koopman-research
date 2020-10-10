function res=boundary_draw(basis,m,n,yt)
load([basis,'_boundary_norepeat_x0.mat'],'X');
for i=1:length(X)
    Y{i}=sort(X{i});
end
hold on
res=drawBoundary(Y,m,n,yt);
%title('Boundary of Logistic Map (x=0)')
end

function res=drawBoundary(X,m,n,yt)
size=10;rate=0.87;
%s=flipud(autumn(n-m+1));
s=jet(n-m+1);
res=[];
for j=m:n
    x=X{j};
    res=[res,x];
    plot(x,yt,'bo','MarkerSize',size*rate^3,'MarkerFaceColor',s(j-m+1,:));%,'Color',s(j-m+1,:)
end
end

