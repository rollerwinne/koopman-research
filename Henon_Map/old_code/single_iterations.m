n=1000;
x =0.3535;
y =0.6869;
k=1.5;
f=@(x,y)mod((x+y+k/(2*pi)*sin(2*pi*x)),1);
g=@(x,y)mod((y+k/(2*pi)*sin(2*pi*x)),1);
X=[x];
Y=[y];
for i=1:n
    x0=f(x,y);
    y0=g(x,y);
    X=[X,x0];
    Y=[Y,y0];
end
Z=[X,Y];
x_n=Z(:,n);
y_n=Z(:,2*n);
F=@(x)(x_n-x);
G=@(y)(y_n-y);
% F=0;
% G=0;
syms xx;
syms yy;
F1 =@(x)diff(F,xx,n);
G1 =@(y)diff(G,yy,n);
x_n1 =F1(x);
y_n1 =G1(y);
% xx=x;
% yy=y;
x=x-x_n/x_n1;
y=y-y_n/y_n1;
for i=1:n
    x0=f(x,y);
    y0=g(x,y);
    X=[X,x0];
    Y=[Y,y0];
end
for j=1:9  
    subplot(3,3,j)
    hold on
    % plot(X,Y,'*')
    scatter3(X,Y,ones(1,length(X)),'filled');%5,zeros(1,length(X)))
%     plot3(X,Y,ones(1,length(X)))
    colorbar
    colormap(jet)
    view(0,90)
    axis equal
    axis([0 1 0 1])
end