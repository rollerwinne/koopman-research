clear all;clc;close all
n=20;q=1;a=1.4;b=0.3;
[x0,y0]=meshgrid(linspace(-1.5,1.5,n));
x=x0(:)';  
y=y0(:)';
f=@(x,y)y+1-a*x.*x;
g=@(x,y)b*x;
dfx=@(x,y)-2*a*x;
dfy=@(x,y)1;
dgx=@(x,y)b;
dgy=@(x,y)0;
dfgxy=@(x,y)[-2*a*x,1;b,0];
for i=1:length(x)
    D=dfgxy(x(i),y(i));
    [F,L]=eig(D);
    L=diag(L);
    if L(1)>=0
        color1{i}='blue';
    else
        color1{i}='red';
    end
    if L(2)>=0
        color2{i}='blue';
    else
        color2{i}='red';
    end
    DX1(i)=F(1,1);
    DY1(i)=F(2,1);
    DX2(i)=F(1,2);
    DY2(i)=F(2,2);
end
hold on 
for i=1:length(x)
    quiver(x(i),y(i),DX1(i),DY1(i),'color',color1{i},'AutoScaleFactor',0.1);
    quiver(x(i),y(i),-DX1(i),-DY1(i),'color',color1{i},'AutoScaleFactor',0.1);
    quiver(x(i),y(i),DX2(i),DY2(i),'color',color2{i},'AutoScaleFactor',0.1);
    quiver(x(i),y(i),-DX2(i),-DY2(i),'color',color2{i},'AutoScaleFactor',0.1);
end
colormap hsv
axis equal