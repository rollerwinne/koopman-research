%% 得到初始值
global n;
global k;
n=2;
x =0.4949;
y =0.5051;
k=1.2;
f=@(x,y)mod((x+y+k/(2*pi)*sin(2*pi*x)),1);
g=@(x,y)mod((y+k/(2*pi)*sin(2*pi*x)),1);
% X=[x];
% Y=[y];
X=zeros(n,1);
Y=zeros(n,1);
for i=1:n
    x=f(x,y);
    y=g(x,y);
    X(i,:)=[x];
    Y(i,:)=[y];
end
%% 求delta
Z=zeros(2*n,1);
J=eye(2*n);
Z=project4_period_F(X,Y,f,g);
J=project4_period_FD(X,Y);
A=J\Z;
A1=sum(A.*A);
delta=sqrt(A1);

if(delta>exp(-6))
    X=X+Z(1:2:end);
    Y=Y+Z(2:2:end);
    Z=project4_period_F(X,Y,f,g);
    J=project4_period_FD(X,Y);
    A=J\Z;
    A1=sum(A.*A);
    delta=sqrt(A1);
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