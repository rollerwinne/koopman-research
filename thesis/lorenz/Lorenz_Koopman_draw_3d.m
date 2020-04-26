clc;close all
clearvars -Except F D K x_function_num
dj=5;%m=4;
choose='complex';
bin=0; % 0:不进行二值化处理 1:二值化处理
tspan=[0,100,0.01];
x0=[-1,3,4];
[t,x]=Lorenz_Koopman_x(tspan,x0);
x_k=x(1:end-1,:);
x_l=x(2:end,:);
%[F,D,K,x_function_num]=Lorenz_Koopman_U(x_k,x_l,dj);
%save('./data/Lorenz_Matrix_data_FD_t0_100_0.01dj5.mat','F','D','K','x_function_num');
%load('./data/Lorenz_Matrix_data_FD_t0_100_0.01dj5.mat');
if bin==1
    F(F>=0)=1;F(F<0)=-1; % binarization
end
if strcmp(choose,'real')==1
    h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.9964 & abs(D)<1.005 & imag(D)>-1e-6 ); % find complex eigenvalues
end
[~,idx]=sort(abs(D(h)),'descend');
h=h(idx);
for i=1:min(9,length(h))
    %figure(i)
    tightsub(3,3,i,0.8);
    fullscreen;
    E_start=1;
    for j=4
        %subplot(2,2,j)
        X=x_k(E_start:end,1);
        Y=x_k(E_start:end,2);
        Z=x_k(E_start:end,3);
        E_temp=K*F(:,h(i));
        if j==1
            E=real(E_temp);
            hh=scatter3(X,Y,Z,3,E,'filled');
            zlabel('real')
        elseif j==2
            E=imag(E_temp);
            hh=scatter3(X,Y,Z,3,E,'filled');
            zlabel('imaginary')
        elseif j==3
            E=abs(E_temp);
            hh=scatter3(X,Y,Z,3,E,'filled');
            zlabel('abs')
        elseif j==4
            E=angle(E_temp);
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
    
    str2=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
    title(str2)
    %suptitle({str1;str2});
    %str=['./temp/Lorenz_eigenfunctions3d_real_figure' num2str(i) '.fig'];
    str=['./temp/Lorenz_eigenfunctions3d01_',choose,'_figure' num2str(i) '.fig'];
    %saveas(hh,str);
    %str=['./temp/Lorenz_eigenfunctions3d_real_figure' num2str(i) '.png'];
    str=['./temp/Lorenz_eigenfunctions3d01_',choose,'_figure' num2str(i) '.png'];
    %saveas(hh,str);
end
str1=[' (n=',num2str(length(x_k(1:end,1))),',m=',num2str(x_function_num(1)),'×',num2str(x_function_num(2)),'×',num2str(x_function_num(3)),...
    '=',num2str(x_function_num(1)*x_function_num(2)*x_function_num(3)),',d_{j}=',num2str(dj),')'];
str='Eigenfunctions of Lorenz System with Gauss Basis';
suptitle([str,str1])
%saveas(gcf,'./temp/Lorenz_eigen_Gauss_leftU_complexL_n10000m960.png');
