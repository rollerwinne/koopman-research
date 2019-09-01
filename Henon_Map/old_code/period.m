%% 得到初始值
global t;
global k;
t=4;
x =0.3838;
y =0.4848;
k=0.02;
f=@(x,y)mod((x+y+k/(2*pi)*sin(2*pi*x)),1);
g=@(x,y)mod((y+k/(2*pi)*sin(2*pi*x)),1);
% X=[x];
% Y=[y];
X=zeros(t,1);
Y=zeros(t,1);
for i=1:t
    x=f(x,y);
    y=g(x,y);
    X(i,:)=[x];
    Y(i,:)=[y];
end
%% 求delta
Z=project4_period_F(X,Y,f,g);
J=project4_period_FD(X);
A=J\(-Z);
A1=sum(A.*A);
delta=sqrt(A1);
count=0;
while count<80
    if (delta>1e-6)
        X=X+0.2*A(1:2:end);
        Y=Y+0.2*A(2:2:end);
        Z=project4_period_F(X,Y,f,g);
        J=project4_period_FD(X);
        A=J\(-Z);
        A1=sum(A.*A);
        delta=sqrt(A1);
        count=count+1;
    else
        break;
    end
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