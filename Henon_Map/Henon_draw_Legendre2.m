%% Initialization
clc;%close all;
clearvars -except F D x_function_length;
tic;timestart=char(datetime('now'));
disp('The running program is from ZC. 么么哒')
%% Parameter settings
n=100;q=1;a=1.4;b=0.3;
power=50;
[x,y]=Henon_x(n,a,b,q);
x_k=[x(1,:);y(1,:)]';
x_l=[x(2,:);y(2,:)]';
xi=abs(x_l(:,1))>1.5 | abs(x_l(:,2))>1.5;
x_k(xi,:)=[];x_l(xi,:)=[];
%% Caculate Eigenfunction
[F,D,x_function_length]=Henon_U_Legendre(x_k,x_l,power);
%save('.\data\Henon_Matrix_data_Legendre_FD_n100power50.mat','F','D','x_function_length');
%load('.\data\Henon_Matrix_data_Legendre_FD_n100power50.mat')
%% Data processing
choose='complex';
if strcmp(choose,'real')==1
    %h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
    h=find(real(D)>0& abs(D)>0.2 & abs(D)<1.5 & abs(imag(D))<1e-6);
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.8 & abs(D)<1.2 & imag(D)>-1e-6 ); % find complex eigenvalues
end
%% Draw Eigenfunctions
figure_num=10;
attachments=[];
for i=1:min(figure_num,length(h))
    figure(i)
    set(gcf,'outerposition',get(0,'screensize'));
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    N=log(b)/log(d_abs);
    A=d_angle/N;
    T=round(360/A);
    err=abs((A*T-360))/360*100;
    str0=['n=',num2str(n),'; power=',num2str(power),'; m=',num2str(x_function_length)];
    str1=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
    str2=['log' '_{' num2str(b) '}(' num2str(d_abs) ')=' num2str(N) '; ' num2str(d_angle) '°/' num2str(N) '=' num2str(A) '°'];
    str3=['T=' num2str(360/A) '≈' num2str(T) '; err=' num2str(err) '%'];
    for j=1:4
        subplot(2,2,j)
        X=real(x_k(:,1));
        Y=real(x_k(:,2));
        if j==1
            Z=real(F(:,h(i)));
        elseif j==2
            Z=imag(F(:,h(i)));
        elseif j==3
            Z=abs(F(:,h(i)));
        elseif j==4
            Z=angle(F(:,h(i)))/pi*180;
        end
        %hh=surf(X,Y,Z);
        hh=scatter3(X(:),Y(:),zeros(length(X(:)),1),7,Z(:),'filled');%scatter3(X,Y,Z,S,C), 前三个参数是坐标，S是散点的size,C是颜色参数
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
    end
    suptitle({str0;[str2,'; ',str3];str1})
    str=['.\temp\Henon_eigenfunctions_Legendre_',num2str(choose),'_figure' num2str(i)];
    saveas(hh,[str,'.fig']);
    saveas(hh,[str,'.png']);
    attachments{i}=[str,'.png'];
end
%% send an E-mail to me
%[subject,content]=qqmail2me(timestart,mfilename('fullpath'),attachments); %程序开始时间、文件名、附件
