function Henon_T_linear_draw(T,Choose,yt)
a=1.4;b=0.3;
f=@(x,y)y+1-a*x.*x;
g=@(x,y)b*x;
period=load('./data/Henon_period_orbrits.mat'); % 周期轨道数据载入
P=period.P;
for choose=Choose
    datastr=['./data/Henon_T_linear_' num2str(T) '_choose' num2str(choose) '.mat'];
    
    if exist(datastr)
        load(datastr);
    else
        [x1,y1,DX1,DX2,DY1,DY2,color1,color2]=T_linear_computing(f,g,P,T,choose);
        save(datastr,'T','x1','y1','DX1','DX2','DY1','DY2','color1','color2');
    end
    hold on
    scale=0.25;
    for i=1:length(x1)
        quiver3(x1(i),y1(i),yt,DX1(i),DY1(i),0,'color',color1{i},'AutoScaleFactor',scale);
        quiver3(x1(i),y1(i),yt,-DX1(i),-DY1(i),0,'color',color1{i},'AutoScaleFactor',scale);
        quiver3(x1(i),y1(i),yt,DX2(i),DY2(i),0,'color',color2{i},'AutoScaleFactor',scale);
        quiver3(x1(i),y1(i),yt,-DX2(i),-DY2(i),0,'color',color2{i},'AutoScaleFactor',scale);
    end
    %     colormap hsv
    axis equal
    %title(['p=' num2str(T)])
    %plot3(x1,y1,0.03*ones(1,length(x1)),'o','color','black','Markersize',4,'MarkerFaceColor','black');
    %drawnow
    % end
end
end

function [x1,y1,DX1,DX2,DY1,DY2,color1,color2]=T_linear_computing(f,g,P,T,choose)
x1=P{T}(choose,mod(1:end,2)==1);
y1=P{T}(choose,mod(1:end,2)==0);
syms x y
for i=1:T
    x_temp=f(x,y);
    y_temp=g(x,y);
    x=x_temp;
    y=y_temp;
end
X=x;
Y=y;
syms x y
J=[diff(X,x),diff(X,y);diff(Y,x),diff(Y,y)];

for i=1:length(x1)
    D=double(subs(J,[x,y],[x1(i),y1(i)]));
    [F,L]=eig(D);
    F=real(F);
    L=diag(L);
    if abs(L(1))<1
        color1{i}='blue';
    else
        color1{i}='red';
    end
    if abs(L(2))<1
        color2{i}='blue';
    else
        color2{i}='red';
    end
    DX1(i)=F(1,1);
    DY1(i)=F(2,1);
    DX2(i)=F(1,2);
    DY2(i)=F(2,2);
end
end