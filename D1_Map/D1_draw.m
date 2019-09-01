%% Initialization
clc;close all;
clear;%vars -except F D;
tic;timestart=char(datetime('now'));
disp('The running program is from ZC. 么么哒')

n=1000;q=1;m=100;md=95;

a=0.3*1i;

x=D1_x(n,a,q);
x_k=x(1,:);
x_l=x(end,:);
%% Caculate Eigenfunction
[F,D] = D1_U(x_k,x_l,m,md);
%[F,D] = D2_U_Fourier(x_k,x_l,y_k,y_l,m);
%save('.\data\Henon_Matrix_data_FD_n100m50md45a1.4b0.3.mat','F','D')
%load('.\data\Henon_Matrix_data_FD_n100m50md45a1.4b0.3.mat')
%% Data processing
choose='complex';
if strcmp(choose,'real')==1
    %h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
    h=find(real(D)>0& abs(D)>0.5 & abs(D)<1.5 & abs(imag(D))<1e-6);
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0 & abs(D)<1.8 & imag(D)>-1e-6 ); % find complex eigenvalues
end
%% Draw Eigenfunctions
figure_num=6;
for i=1:min(figure_num,length(h))
    figure(i)
    set(gcf,'outerposition',get(0,'screensize'));
    for j=1:4
        subplot(2,2,j)
        X=real(x_k);
        if j==1
            Y=real(F(:,h(i)));
        elseif j==2
            Y=imag(F(:,h(i)));
        elseif j==3
            Y=abs(F(:,h(i)));
        elseif j==4
            Y=angle(F(:,h(i)))/pi*180;
        end
        hh=stem(X,Y,'.');
        %hh=scatter3(X(:),Y(:),zeros(length(X(:)),1),3,Z(:));%scatter3(X,Y,Z,S,C), 前三个参数是坐标，S是散点的size,C是颜色参数
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
        xlim([-2,2])
    end
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    
    str0=['a=',complex2str(a),' n=',num2str(n),'; m=',num2str(m),'; dj=4/',num2str(md)];
    str1=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
    suptitle({str0;str1})
    
    str5=['a_',complex2str(a)];
    str6=['D1_eigenfunctions_',num2str(choose),'_',str5,'_figure' num2str(i)];
    str7='.\fig\';
    warning off
    mkdir(str7,str5);
    warning on
    str8=[str7,str5,'\',str6];
    saveas(hh,[str8,'.fig']);
    saveas(hh,[str8,'.png']);
end
