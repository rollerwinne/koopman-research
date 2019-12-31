clear
a=1.4;b=0.3;n=2000;nn=50;
m=20;md=m;

[x,y]=meshgrid(linspace(-1.5,1.5,nn));
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
f_inv=@(x,y)y/b;
g_inv=@(x,y)x-1+a/b/b*y.*y;
x_fix=((b-1)+((b-1).^2+4*a).^(1/2))*((2*a).^-1);
y_fix=b*x_fix;
A=[-2*a*x_fix,1;b,0];
[V,~]=eig(A);                     %计算本征向量
k=V(2,1)/V(1,1);
x_step=linspace(-0.001+2e-4,0.002,n);
x0=x_fix+x_step;
y0=k*(x0-x_fix)+y_fix;

times=12;
for i=1:times
    x_temp=f(x0,y0);
    y_temp=g(x0,y0);
    x0=x_temp;
    y0=y_temp;
    % figure
    % plot(x0,y0,'b.')
end
for i=1:30
    x=f(x,y);
    y=g(x,y);
end
x(isinf(x))=[];
y(isinf(y))=[];
[~,idx]=sort(y0);
x0=x0(idx);
y0=y0(idx);

x_k=x(:);
y_k=y(:);
% x_k=x0(idx);
% y_k=y0(idx);

x_l=f(x_k,y_k);
y_l=g(x_k,y_k);
%% Data Reading
Attr=load('./data/Henon_attractors_data_xy.mat'); % 吸引子数据载入
x_attr=Attr.x;y_attr=Attr.y;
Peri=load('./data/Henon_period_orbrits_P_1_0.3_-1_-0.3.mat');
P=Peri.P;
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
    h=find(abs(D)>0.1 & abs(D)<1.05 & imag(D)>-1e-6 ); % find complex eigenvalues
    h=1:10;
end
%% Draw Eigenfunctions
figure_num=10;
attachments=[];
for i=1:min(figure_num,length(h))
    figure(i)
    set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
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
        X=real(reshape(x_k,1,length(x_k)));
        Y=real(reshape(y_k,1,length(x_k)));
        if j==1
            Z=real(reshape(F(:,h(i)),1,length(x_k)));
        elseif j==2
            Z=imag(reshape(F(:,h(i)),1,length(x_k)));
        elseif j==3
            Z=abs(reshape(F(:,h(i)),1,length(x_k)));
        elseif j==4
            Z=angle(reshape(F(:,h(i)),1,length(x_k)))/pi*180;
        end
        res=find_similiar(x0,y0,X,Y,Z);
        X=res(:,1);
        Y=res(:,2);
        Z=res(:,3);
        %hh=surf(X,Y,Z);
        %hh=scatter3(X(:),Y(:),zeros(length(X(:)),1),3,Z(:));%scatter3(X,Y,Z,S,C), 前三个参数是坐标，S是散点的size,C是颜色参数
        scatter3(X(:),Y(:),Z(:),3,Z(:));
        hold on
        z_min=min(Z(:));z_max=max(Z(:));
        z_min1=z_min+0.01*(z_max-z_min);
%         scatter3(x_attr(:),y_attr(:),z_min*ones(length(x_attr(:)),1),3,min(Z(:))*ones(length(x_attr(:)),1));%画吸引子
%         plot3([Attr.x1;Attr.x2],[Attr.y1;Attr.y2],z_min1*ones(2,1),'r*');
        
        %B=Henon_Boundary(xy_bound,5,5);
        %plot3(B(:,1),B(:,2),z_min1*ones(length(B(:,1)),1),'ro');
        load('boundary.mat','boundary')
        hh=plot3(boundary(1,:),boundary(2,:),z_min1*ones(length(boundary(1,:)),1),'ro');
       
        view(-15,60)
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        if j==1
            ylabel('real')
        elseif j==2
            ylabel('imaginary')
        elseif j==3
            ylabel('abs')
        elseif j==4
            ylabel('angle')
        end
        shading interp
        colorbar
        colormap(jet)
    end
    suptitle({str0;[str2,'; ',str3];str1})
    % str=['./temp3/Henon_eigen_attr_',num2str(choose),'_figure' num2str(i)];
    % saveas(hh,[str,'.fig']);
    % saveas(hh,[str,'.png']);
    % attachments{i}=[str,'.png'];
end

function res=find_similiar(x0,y0,X,Y,Z)
res=zeros(length(X),3);
num=1;
for i=1:length(X)
    for j=1:length(x0)
        if ((X(i)-x0(j))^2+(Y(i)-y0(j))^2)^0.5<0.05
            res(num,:)=[X(i),Y(i),Z(i)];
            num=num+1;
            break;
        end
    end
end
res(num:end,:)=[];
end