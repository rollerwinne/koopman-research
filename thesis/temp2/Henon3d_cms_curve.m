clear;close all;clc
%% flow
%选取fixed point旁边的一组小线段

A=-1.86;B=0.72;C=0.03;
x(:,1)=0.1; y(:,1)=0;z(:,1)=0;
for i=1:10000
   x(:,i+1)=y(:,i);
   y(:,i+1)=z(:,i);
   z(:,i+1)=-1.45.*z(:,i).^2+0.515.*y(:,i).*z(:,i)-y(:,i).^2+B.*x(:,i)+C.*y(:,i)+A.*z(:,i);                                  %-y(:,i).^2+a.*x(:,i)+b.*z(:,i)+c;
end
figure(1)
plot3(x,y,z,'b.')
xlabel('x_n')
ylabel('y_n')
zlabel('z_n')

%% 寻找不动点
global F G H F1 G1 H1  
F=@(x,y,z)y;
G=@(x,y,z)z;
H=@(x,y,z)B.*x+C.*y+A.*z-y.^2+0.515.*y.*z-1.45.*z.^2;
F1=@(x,y,z)1/B.*(z-(C.*x+A.*y-1.45.*y.^2+0.515.*x.*y-x.^2));
G1=@(x,y,z)x;
H1=@(x,y,z)y;
syms symx;
equ=(0.515-1-1.45).*symx.^2+(B+C+A-1).*symx;
[ansx]=solve(equ,symx);
x_fix=double(ansx);
y_fix=x_fix;
z_fix=x_fix;
figure(2)
plot3(x_fix,y_fix,z_fix,'r*')
hold on
plot3(x,y,z,'b.')

col=2;
A_matrix=[0,1,0;
    0,0,1;
    B,C-2.*y_fix(col,:)+0.515.*z_fix(col,:),A+0.515.*y_fix(col,:)-2.9.*z_fix(col,:)];
[V,D]=eig(A_matrix); %计算本征向量  
x_step=-0.001:0.0000005:0.001;
x0=x_fix(col,:)+x_step;
y0=V(2,1).*(x0-x_fix(col,:))/V(1,1)+y_fix(col,:);%计算本征向量所对应直线，得到线段上一组点
z0=V(3,1).*(x0-x_fix(col,:))/V(1,1)+z_fix(col,:);
figure(3)
plot3(x0,y0,z0,'r*')
hold on 
plot3(x,y,z,'b.')
xlabel('x_n')
ylabel('y_n')
zlabel('z_n')
title('Original data')

%小线段不断演化得到更大轨迹线
bag{1}=[x0;y0;z0]
iteration=40;
for i=1:iteration
x0_evo=F(x0,y0,z0);
y0_evo=G(x0,y0,z0);
z0_evo=H(x0,y0,z0);
bag{i+1}=[x0_evo;y0_evo;z0_evo];
x0=x0_evo;
y0=y0_evo;
z0=z0_evo;
%figure(3+i)
%plot3(x0,y0,z0,'b-')
%xlabel('x_n')
%ylabel('y_n')
%zlabel('z_n')
%aa=i;
%title(['N=',num2str(aa)])
%axis([-1.5,0.5,-1.5,0.5,-1.5,0.5]) 
end
%save('bag')
figure(4)
plot3(bag{17}(1,:),bag{17}(2,:),bag{17}(3,:),'b-')
xlabel('x_n')
ylabel('y_n')
zlabel('z_n')
axis([-1.5,0.5,-1.5,0.5,-1.5,0.5]) 
x_boundary1=bag{17}(1,:);
y_boundary1=bag{17}(2,:);
z_boundary1=bag{17}(3,:);
l_boundary1=length(x_boundary1);
sp1 = spline(1:l_boundary1,[x_boundary1;y_boundary1;z_boundary1],1:.1:l_boundary1); 
clear Curve n_curve dif_p mp1 n_mp1 nb n_sym1 n_useful 
Curve=curve3(sp1(1,:),sp1(2,:),sp1(3,:));    %调用曲率函数，计算曲率,找到曲率极大值
t_step=2/(length(Curve)-1);
ti=-1:t_step:1;
dcurve=gradient(Curve)/t_step;
[~,n_max]=find(diff(sign(diff(Curve)))==-2);

figure
%plot3(x_boundary1,y_boundary1,z_boundary1,'r-',sp1(1,:),sp1(2,:),sp1(3,:),'linewidth',1.5)
xlabel('x_n')
ylabel('y_n')
zlabel('z_n')
axis([-1.5,0.5,-1.5,0.5,-1.5,0.5])
plot3(sp1(1,:),sp1(2,:),sp1(3,:),'b-')
hold on
plot3(sp1(1,n_max),sp1(2,n_max),sp1(3,n_max),'r*')
hold on
plot3(x_fix(col,:),y_fix(col,:),z_fix(col,:),'k+')
%plot3(sp1(1,n_max(:,1)),sp1(2,n_max(:,1)),sp1(3,n_max(:,1)),'k*')
hold on
n_kl(1,:)=Curve(:,n_max);
%n_kl(2,:)=n_max;

%xx_test=sp1(1,n_max(:,4));
%yy_test=sp1(2,n_max(:,4));
%zz_test=sp1(3,n_max(:,4));
%xx_evo=F(xx_test,yy_test,zz_test);
%yy_evo=G(xx_test,yy_test,zz_test);
%zz_evo=H(xx_test,yy_test,zz_test);
%plot3(xx_evo,yy_evo,zz_evo,'y*')

%得到base line
p_data0=sp1(:,n_max(:,2):n_max(:,end-1));    %得第一条轨迹线，这里如何选取边界点需极为谨慎！！！
x_boundary2=p_data0(1,:);
y_boundary2=p_data0(2,:);
z_boundary2=p_data0(3,:);
l_boundary2=length(x_boundary2);
step0=l_boundary2/40000;
sp2 = spline(1:l_boundary2,[x_boundary2;y_boundary2;z_boundary2],1:step0:l_boundary2); 
global p_data1
p_data1=sp2;
N=length(p_data1);
figure
plot3(p_data1(1,:),p_data1(2,:),p_data1(3,:),'b-')
xlabel('x_n')
ylabel('y_n')
zlabel('z_n')
axis([-1.5,0.5,-1.5,0.5,-1.5,0.5]) %得到初步边界轨迹

%baseline不断演化，得到后续轨迹线
x0=sp2(1,:);
y0=sp2(2,:);
z0=sp2(3,:);
figure
for i=1:3
x0_evo=F(x0,y0,z0);
y0_evo=G(x0,y0,z0);
z0_evo=H(x0,y0,z0);
bag{i+1}=[x0_evo;y0_evo;z0_evo];
x0=x0_evo;
y0=y0_evo;
z0=z0_evo;
plot3(x0,y0,z0,'b-')
hold on
xlabel('x_n')
ylabel('y_n')
zlabel('z_n')
aa=i;
title(['N=',num2str(aa)])
axis([-1.5,0.5,-1.5,0.5,-1.5,0.5]) 
end






