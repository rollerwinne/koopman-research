%% Initialization
clc;close all;
clearvars -except F D;
tic;timestart=char(datetime('now'));
%disp('The running program is from ZC. 么么哒')
%% Parameter settings
a=1.4;b=0.3;q=50;n=40;
m=20;md=m;
Attr=load('./data/Henon_attractors_data_xy.mat'); % 吸引子数据载入
x_attr=Attr.x;y_attr=Attr.y;
Peri=load('./data/Henon_period_orbrits_P_1_0.3_-1_-0.3.mat');
P=Peri.P;
xy_bound=[1.272933828112852,-0.012403304471461];
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
[x0,y0]=meshgrid(linspace(-1.5,1.5,n));
for i=1:q
    x_temp=f(x0,y0);
    y_temp=g(x0,y0);
    x0=x_temp;
    y0=y_temp;
end
x_k=x0(:)';
y_k=y0(:)';
x_l=f(x_k,y_k);
y_l=g(x_k,y_k);
%% Caculate Eigenfunction
[F,D,U] = Henon_U(x_k,x_l,y_k,y_l,m,md);
%save('.\data\Henon_Matrix_data_on_attractors_FD_n50m50.mat','F','D');
%load('.\data\Henon_Matrix_data_on_attractors_FD_n50m50.mat');
%% Data processing
choose='complex';
if strcmp(choose,'real')==1
    %h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
    h=find(real(D)>0& abs(D)>0.0001 & abs(D)<1.2 & abs(imag(D))<1e-6);
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.5 & abs(D)<1.05 & imag(D)>-1e-6 ); % find complex eigenvalues
end
%% Draw Eigenfunctions
figure_num=10;
attachments=[];
for i=1:min(figure_num,length(h))
    figure(i)
    set(gcf,'outerposition',get(0,'screensize'));
    %ha = tight_subplot(1,2,[0 .03],[.2 .1],[.02 .02]);
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    N=log(b)/log(d_abs);
    A=d_angle/N;
    T=round(360/A);
    err=abs((A*T-360))/360*100;
    str0=['n=',num2str(n),'; m=',num2str(m),'; dj=3/',num2str(md)];
    str1=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
    str2=['log' '_{' num2str(b) '}(' num2str(d_abs) ')=' num2str(N) '; ' num2str(d_angle) '°/' num2str(N) '=' num2str(A) '°'];
    str3=['T=' num2str(360/A) '≈' num2str(T) '; err=' num2str(err) '%'];
    for j=1:4
        %figure
        subplot(2,2,j)
        %axes(ha(j));
        X=real(reshape(x_k,n,n));
        Y=real(reshape(y_k,n,n));
        if j==1
            Z=real(reshape(F(:,h(i)),n,n));
        elseif j==2
            Z=imag(reshape(F(:,h(i)),n,n));
        elseif j==3
            Z=abs(reshape(F(:,h(i)),n,n));
        elseif j==4
            Z=angle(reshape(F(:,h(i)),n,n))/pi*180;
        end
        %hh=surf(X,Y,Z);
        %hh=scatter3(X(:),Y(:),zeros(length(X(:)),1),3,Z(:));%scatter3(X,Y,Z,S,C), 前三个参数是坐标，S是散点的size,C是颜色参数
        scatter3(X(:),Y(:),Z(:),3,Z(:));
        hold on
        z_min=min(Z(:));z_max=max(Z(:));
        z_min1=z_min+0.01*(z_max-z_min);
        scatter3(x_attr(:),y_attr(:),z_min*ones(length(x_attr(:)),1),3,min(Z(:))*ones(length(x_attr(:)),1));%画吸引子
        plot3([Attr.x1;Attr.x2],[Attr.y1;Attr.y2],z_min1*ones(2,1),'r*');
        
        %B=Henon_Boundary(xy_bound,5,5);
        %plot3(B(:,1),B(:,2),z_min1*ones(length(B(:,1)),1),'ro');
        load('boundary.mat','boundary')
        hh=plot3(boundary(1,:),boundary(2,:),z_min1*ones(length(boundary(1,:)),1),'ro');
        
        view(-15,60)
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        if j==1
            ylabel('real')
            %title({str1;str2;str3})
        elseif j==2
            ylabel('imaginary')
            %title(str2)
        elseif j==3
            ylabel('abs')
            %title(str3)
        elseif j==4
            ylabel('angle')
        end
        shading interp
        colorbar
        colormap(jet)
    end
    suptitle({str0;[str2,'; ',str3];str1})
    str=['./temp2/Henon_eigen_attr_',num2str(choose),'_figure' num2str(i)];
    % saveas(hh,[str,'.fig']);
    saveas(hh,[str,'.png']);
    % attachments{i}=[str,'.png'];
end
%% send an E-mail to me
%[subject,content]=qqmail2me(timestart,mfilename('fullpath'),attachments); %程序开始时间、文件名、附件