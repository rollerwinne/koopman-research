clear;clc;close all;
A=-1.86;B=0.72;C=0.03;
f=@(x,y,z)y;g=@(x,y,z)z;h=@(x,y,z)-1.45.*z.^2+0.515.*y.*z-y.^2+B.*x+C.*y+A.*z;
fun=@(x)[x(2);x(3);B*x(1)+C*x(2)+A*x(3)-1.45*x(3).^2+0.515.*x(2).*x(3)-x(2).^2];
fun2=@(x)[f(x(1),x(2),x(3));g(x(1),x(2),x(3));h(x(1),x(2),x(3))];
f_inv=@(x,y,z)1/B.*(z-(C.*x+A.*y-1.45.*y.^2+0.515.*x.*y-x.^2));g_inv=@(x,y,z)x;h_inv=@(x,y,z)y;
fun_inv=@(x)[1/B.*(x(3)-(C.*x(1)+A.*x(2)-1.45.*x(2).^2+0.515.*x(1).*x(2)-x(1).^2));x(1);x(2)];
fun_inv2=@(x)[f_inv(x(1),x(2),x(3));g_inv(x(1),x(2),x(3));h_inv(x(1),x(2),x(3))];

syms x;
equ=(0.515-1-1.45).*x.^2+(B+C+A-1).*x;
sx=solve(equ,x);
x_fix=double(sx(2));
y_fix=x_fix;
z_fix=x_fix;
J=[0,1,0;
    0,0,1;
    B,0.515.*z_fix-2.*y_fix+C,-1.45*2.*z_fix+0.515.*y_fix+A];

x0=[0.1,0,0];
iter=10000;
x=x0;
for i=1:iter
    x=fun(x);
    X(:,i)=x;
end
%plot3(X(1,:),X(2,:),X(3,:),'b.')

[V,D]=eig(J); %计算本征向量
x_step=-0.01:0.0000001:0.01;
x1=x_fix+x_step';
y1=V(2,1).*(x1-x_fix)/V(1,1)+y_fix;%计算本征向量所对应直线，得到线段上一组点
z1=V(3,1).*(x1-x_fix)/V(1,1)+z_fix;

for i=1:29
    Y(:,:,i)=[x1,y1,z1];%1点的数量2三维3第几次演化
    x_temp=f(x1,y1,z1);
    y_temp=g(x1,y1,z1);
    z_temp=h(x1,y1,z1);
    x1=x_temp;
    y1=y_temp;
    z1=z_temp;
    %     plot3(Y(:,1,i),Y(:,2,i),Y(:,3,i),'b');
    %     pause(0.05);
    %     title(i)
    %     drawnow;
end

x2=f(x1,y1,z1);y2=g(x1,y1,z1);z2=h(x1,y1,z1);
x3=f(x2,y2,z2);y3=g(x2,y2,z2);z3=h(x2,y2,z2);

dim=3;
m=2;n=length(x1);
K=zeros(n*dim,m);L=K;
K(1:dim:end,1)=x1;K(2:dim:end,1)=y1;K(3:dim:end,1)=z1;
L(1:dim:end,1)=x2;L(2:dim:end,1)=y2;L(3:dim:end,1)=z2;
K(1:dim:end,2)=x2;K(2:dim:end,2)=y2;K(3:dim:end,2)=z2;
L(1:dim:end,2)=x3;L(2:dim:end,2)=y3;L(3:dim:end,2)=z3;
U=pinv(K)*L;
[F,D]=eig(U);
D=diag(D);
%h=1:length(D);
hh=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6);
if length(D)<=9
    hh=1:length(D);
end
A=K*F(:,hh(1));
A_real=real(A);
A_abs=abs(A);

x_str=['x','y','z'];
figure;
set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
for i=1:3
    subplot(2,dim,i)
    idx=findp(A_real(i:dim:end));
    hold on;
    %scatter3(x1,y1,z1,3,A_real(i:dim:end));
    plot3(x1,y1,z1);
    plot3(x1(idx),y1(idx),z1(idx),'o','Color','red','MarkerFaceColor','red');
    xlabel('x');ylabel('y');zlabel('z');
    xlabel([x_str(i),'-real']);
    Henon_cms_draw('black',f_inv,g_inv,h_inv);
    graph_control_3D;
    
    subplot(2,dim,i+dim)
    idx=findp(A_abs(i:dim:end));
    %scatter3(x1,y1,z1,3,A_abs(i:dim:end));
    hold on;
    plot3(x1,y1,z1);
    plot3(x1(idx),y1(idx),z1(idx),'o','Color','red','MarkerFaceColor','red');
    xlabel([x_str(i),'-abs']);
    xlabel('x');ylabel('y');zlabel('z');
    Henon_cms_draw('black',f_inv,g_inv,h_inv);
    graph_control_3D;
end

function graph_control_3D
view([0,0])
colorbar
colormap(jet)
end

function idx=findp(A)
a=max(A)-min(A);
[~,index1]=findpeaks(A,'minpeakheight',a/20,'minpeakdistance',length(A)/20);
[~,index2]=findpeaks(-A,'minpeakheight',a/20,'minpeakdistance',length(A)/20);
idx=[index1;index2];
end

