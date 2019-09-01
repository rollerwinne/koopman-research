%% Initialization
clc;%close all;
clearvars -except F D;
tic;timestart=char(datetime('now'));
disp('The running program is from ZC. 么么哒')
%% Parameter settings
a=1.4;b=0.3;q=50;n=40;
m=8;md=m;
%load('.\data\Henon_attractors_data_xy.mat'); % 吸引子数据载入
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
[x,y]=meshgrid(linspace(-1.5,1.5,n));
for i=1:q
    x_temp=f(x,y);
    y_temp=g(x,y);
    x=x_temp;
    y=y_temp;
end
x_k=x(:)';
y_k=y(:)';
x_l=f(x_k,y_k);
y_l=g(x_k,y_k);
%% Caculate Eigenfunction
[F,D,U] = Henon_U(x_k,x_l,y_k,y_l,m,md);
%save('.\data\Henon_Matrix_data_on_attractors_FD_n50m50.mat','F','D');
%load('.\data\Henon_Matrix_data_on_attractors_FD_n50m50.mat');
%% Data processing
choose='real';
if strcmp(choose,'real')==1
    %h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
    h=find(real(D)>0& abs(D)>0.0001 & abs(D)<1.2 & abs(imag(D))<1e-6);
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.95 & abs(D)<1.05 & imag(D)>-1e-6 ); % find complex eigenvalues
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
    str0=['n=',num2str(m),'; m=',num2str(m),'; dj=3/',num2str(md)];
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
        hh=scatter3(X(:),Y(:),zeros(length(X(:)),1),3,Z(:));%scatter3(X,Y,Z,S,C), 前三个参数是坐标，S是散点的size,C是颜色参数
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
        view(0,90)
        axis([-1.5 1.5 -1.5 1.5])
        axis equal
        %title(str1)
        %xlabel(str2)
    end
    suptitle({str0;[str2,'; ',str3];str1})
    str=['.\temp\Henon_eigenfunctions_on_attractors_',num2str(choose),'_figure' num2str(i)];
    saveas(hh,[str,'.fig']);
    saveas(hh,[str,'.png']);
    attachments{i}=[str,'.png'];
end
%% send an E-mail to me
%[subject,content]=qqmail2me(timestart,mfilename('fullpath'),attachments); %程序开始时间、文件名、附件