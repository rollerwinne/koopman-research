%% Initialization
clc;%close all;
clearvars -except F D x_function_length;
tic;timestart=char(datetime('now'));
disp('The running program is from ZC. 么么哒')
%% parameter settings
basis_function='legendre';
power=20;
tspan=[0,100,0.01];
x0=[-1,3,4];
[t,x]=Lorenz_Koopman_x(tspan,x0);
x_k=x(1:end-1,:);
x_l=x(2:end,:);
%% Caculate Eigenfunction
[F,D,x_function_length]=Lorenz_Koopman_U_legendre(x_k,x_l,power);
%save('.\data\Lorenz_Matrix_data_legendre_FD_t0_100_0.01power20.mat','F','D','x_function_length');
%load('.\data\Lorenz_Matrix_data_legendre_FD_t0_100_0.01power20.mat');
%% data processing
choose='real'; % real or complex
bin=0; % 0:不进行二值化处理 1:二值化处理
if bin==1
    F(abs(F)>=0)=1;F(abs(F)<0)=-1; % binarization
end
if strcmp(choose,'real')==1
    h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.997 & abs(D)<1.003 & imag(D)>-1e-6 ); % find complex eigenvalues
end
%% draw Eigenfunctions
figure_num=10;
attachments=[];
for i=1:min(figure_num,length(h))
    figure(i)
    set(gcf,'outerposition',get(0,'screensize'));
    E_start=1;
    for j=1:4
        subplot(2,2,j)
        X=x_k(E_start:end,1);
        Y=x_k(E_start:end,2);
        Z=x_k(E_start:end,3);
        if j==1
            E=real(F(E_start:end,h(i)));
            hh=scatter3(X,Y,Z,3,E,'filled');
            zlabel('real')
        elseif j==2
            E=imag(F(E_start:end,h(i)));
            hh=scatter3(X,Y,Z,3,E,'filled');
            zlabel('imaginary')
        elseif j==3
            E=abs(F(E_start:end,h(i)));
            hh=scatter3(X,Y,Z,3,E,'filled');
            zlabel('abs')
        elseif j==4
            E=angle(F(E_start:end,h(i)));
            hh=scatter3(X,Y,Z,3,E,'filled');
            zlabel('angle')
        end
        colorbar
        colormap(jet)
        view([1 -1 0])
        axis equal        
    end
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    str1=['n=',num2str(length(x_k(1:end,1))),'; power=',num2str(power),'; m=',num2str(x_function_length)];
    str2=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
    suptitle({str1;str2});
    str=['.\temp\Lorenz_eigenfunctions_',basis_function,'_',choose,'_figure' num2str(i)];
    saveas(hh,[str,'.fig']);
    saveas(hh,[str,'.png']);
    attachments{i}=[str,'.png'];
end
%% send an E-mail to me
%[subject,content]=qqmail2me(timestart,mfilename('fullpath'),attachments); %程序开始时间、文件名、附件