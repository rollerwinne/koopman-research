clear all;clc;close all
a=1.4;b=0.3;q=1;n=50;
load('.\data\Henon_attractors_data_xy.mat'); % 吸引子数据载入
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
[x,y]=meshgrid(linspace(-1.5,1.5,n));
for i=1:50
    x_temp=f(x,y);
    y_temp=g(x,y);
    x=x_temp;
    y=y_temp;
end
x_k=x(:)';
y_k=y(:)';
x_l=f(x_k,y_k);
y_l=g(x_k,y_k);
for m=20
    for md=15
        [F,D,U] = Henon_U(x_k,x_l,y_k,y_l,m,md);
        %[F,D,U] = project4_U_f(x_k,x_l,y_k,y_l,m);
        %figure('NumberTitle','off','Name',['m=' num2str(m) '  dj=3/' num2str(md)]);
%         h=find(abs(D)<1.5 & imag(D)>-1e-6);
        h=find(real(D)>0&abs(D)<1.2 & abs(imag(D))<1e-6);
%         tic
        for i=1:min(9,length(h))
            subplot(3,3,i)
            %  stem(x0,F(:,h(i)),'.');
            hh=scatter3(x_k,y_k,F(:,h(i))',3,F(:,h(i))');%scatter3(X,Y,Z,S,C), 前三个参数是坐标，S是散点的size,C是颜色参数
            %  caxis([cc(i,1),cc(i,2)]);
            colorbar
            colormap(jet)
            view(0,90)
            xlabel('x');ylabel('y')
            axis([-1.5 1.5 -1.5 1.5])
            axis equal
            title([num2str(abs(D(h(i)))) ' ∠' num2str(angle(D(h(i)))/pi*180) '°'])
         end
    end
end