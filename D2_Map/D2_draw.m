%% Initialization
clc;%close all;
clear;%vars -except F D;
tic;timestart=char(datetime('now'));
disp('The running program is from ZC. 么么哒')
%% Parameter settings
% A=[-1.2,-0.7,0.4,1.6];
% B=[-1.3,-0.6,0.3,1.5];
% [AA,BB]=meshgrid(A,B);
% AA=AA(:);BB=BB(:);
% ABDT=[AA,BB,AA.*BB,AA+BB,(AA+BB).^2-abs(4.*AA.*BB)];

% DD=[0.5,-0.5,0.5,0.5, 1.5,-1.5,	1.5,-1.5, 4 , 1, 1]';
% TT=[2,2,1.2,-1.2,     3,3,        2,2,      4 , 0, 1]';
DD=[0.5,0.5]';
TT=[1.2,-1.2]';
ABDT=[(TT+sqrt(TT.^2-4.*DD))./2,(TT-sqrt(TT.^2-4.*DD))./2,DD,TT,abs(TT).^2-abs(4*DD)];
%%
n=50;q=1;m=20;md=18;
for abi=1:size(ABDT,1)
    close all
    a=ABDT(abi,1);b=ABDT(abi,2);
    det=ABDT(abi,3);trace=ABDT(abi,4);
    
    [x,y]=D2_x(n,a,b,q);
    x_k=x(1,:);
    y_k=y(1,:);%p时刻数据
    x_l=x(end,:);
    y_l=y(end,:);
    %% Caculate Eigenfunction
    [F,D] = D2_U(x_k,x_l,y_k,y_l,m,md);
    %[F,D] = D2_U_Fourier(x_k,x_l,y_k,y_l,m);
    %save('.\data\Henon_Matrix_data_FD_n100m50md45a1.4b0.3.mat','F','D')
    %load('.\data\Henon_Matrix_data_FD_n100m50md45a1.4b0.3.mat')
    %% Data processing
    choose='complex';
    if strcmp(choose,'real')==1
        %h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
        h=find(real(D)>0& abs(D)>0.5 & abs(D)<1.5 & abs(imag(D))<1e-6);
    elseif strcmp(choose,'complex')==1
        h=find(abs(D)>0.8 & abs(D)<1.2 & imag(D)>-1e-6 ); % find complex eigenvalues
    end
    %% Draw Eigenfunctions
    figure_num=6;
    for i=1:min(figure_num,length(h))
        figure(i)
        set(gcf,'outerposition',get(0,'screensize'));
        for j=1:4
            subplot(2,2,j)
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
            hh=surf(X,Y,Z);
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
            shading interp
            colorbar
            colormap(jet)
            view(0,90)
            axis([-2,2,-2,2])
            axis equal
        end
        d_abs=abs(D(h(i)));
        d_angle=angle(D(h(i)))/pi*180;
        N=log(b)/log(d_abs);
        A=d_angle/N;
        T=round(360/A);
        err=abs((A*T-360))/360*100;
        str0=['n=',num2str(n),'; m=',num2str(m),'; dj=2/',num2str(md)];
        str1=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
        str3=['a=',complex2str(a),' b=',complex2str(b),' det=',num2str(det),' trace=',num2str(trace)];
        suptitle({str3;str0;str1})
        
        str5=['a_',complex2str(a),'_b_',complex2str(b),'_det_',num2str(det),'_trace_',num2str(trace)];
        str6=['D2_eigenfunctions_',num2str(choose),'_',str5,'_figure' num2str(i)];
        str7='.\fig\';
        warning off
        mkdir(str7,str5);
        warning on
        str8=[str7,str5,'\',str6];
        saveas(hh,[str8,'.fig']);
        saveas(hh,[str8,'.png']);
    end
    disp(['abi=',num2str(abi)]);
end