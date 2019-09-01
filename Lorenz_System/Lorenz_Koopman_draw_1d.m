clear all;clc;close all
load('.\data\Lorenz_Matrix_data_FD_t0_100_0.01dj5.mat');
dj=5;%m=4;
choose='real';  % 'real' or 'complex'
bin=0;
tspan=[0,100,0.01];
x0=[-1,3,4];
[t,x]=Lorenz_Koopman_x(tspan,x0);
x_k=x(1:end-1,:);
x_l=x(2:end,:);
%[F,D,U,x_function_num]=Lorenz_Koopman_U(x_k,x_l,dj);
%save('.\data\Lorenz_Matrix_data_FD_t0_100_0.02dj5.mat','F','D','x_function_num');
%load('.\data\Lorenz_Matrix_data_FD_t0_100_0.01dj5.mat');
if bin==1
    F(F>=0)=1;F(F<0)=-1; % binarization
end
if strcmp(choose,'real')==1
    h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.95 & abs(D)<1.05 & imag(D)>-1e-6 ); % find complex eigenvalues
end
for i=1:min(4,length(h))
    %subplot(4,1,i)
    figure(i)
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    str2=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
    for j=1:4
        subplot(4,1,j)
        E_start=200;
        if j==1
            E=real(F(:,h(i)));
            hh=stem(t(E_start:end-1),E(E_start:end),'.');
            ylabel('real')
        elseif j==2
            E=imag(F(:,h(i)));
            hh=stem(t(E_start:end-1),E(E_start:end),'.');
            ylabel('imaginary')
        elseif j==3
            E=abs(F(:,h(i)));
            hh=stem(t(E_start:end-1),E(E_start:end),'.');
            ylabel('abs')
        elseif j==4
            E=angle(F(:,h(i)));
            hh=stem(t(E_start:end-1),E(E_start:end),'.');
            ylabel('angle')
        end
    end
    str1=['n=',num2str(length(x_k(1:end,1))),'; m=',num2str(x_function_num(1)),'*',num2str(x_function_num(2)),'*',num2str(x_function_num(3)),...
        '=',num2str(x_function_num(1)*x_function_num(2)*x_function_num(3)),'; d_{j}=',num2str(dj),''];
    suptitle({str1;str2})
    str=['.\temp\Lorenz_eigenfunctions1d_',choose,'_figure' num2str(i) '.fig'];
    %str=['.\temp\Lorenz_eigenfunctions1d_complex_figure' num2str(i) '.fig'];
    saveas(hh,str);
    str=['.\temp\Lorenz_eigenfunctions1d_',choose,'_figure' num2str(i) '.png'];
    %str=['.\temp\Lorenz_eigenfunctions1d_complex_figure' num2str(i) '.png'];
    saveas(hh,str);
end
