%% Initialization
clc;close all;
clearvars -except F D K;
%% Parameter settings
n=100;q=1;a=1.4;b=0.3;
m=50;md=45;
setup.function='Gauss';setup.gauss.m=m;setup.gauss.md=md;
% setup.function='Poly';setup.poly.power=3;
% setup.function='Legendre';setup.legendre.power=3;
% setup.function='Fourier';setup.fourier.m=50;
setup.leftU=false;

%% Caculate Eigenfunction
[x,y]=Henon_x(n,a,b,q);
x_k=x(1,:);y_k=y(1,:);
x_l=x(end,:);y_l=y(end,:);

% attr=load('./data/Henon_attractors_n6245.mat');
% x_k=attr.x(:);y_k=attr.y(:);
% f=@(x,y)y+1-a.*x.*x;
% g=@(x,y)b*x;
% x_l=f(x_k,y_k);y_l=g(x_k,y_k);
% x=[x_k,x_l]';
% y=[y_k,y_l]';
% n=length(x_k);
%[F,D,K] = Henon_basis_FDK(x,y,setup);
%save('./data/Henon_FDK_leftU_n100m50md45a1.4b0.3.mat','F','D','K')
%save('./data/Henon_FDK_n100m50md45a1.4b0.3.mat','F','D','K')
%save('./data/Henon_FDK_Fourier_n100m50a1.4b0.3.mat','F','D','K')
%load('Henon_FDK_n100m50md45a1.4b0.3.mat')
%load('./data/Henon_FDK_leftU_n100m50md45a1.4b0.3.mat','F','D','K');
%% Data processing
choose='complex';
if strcmp(choose,'real')==1
    h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
    %h=find(real(D)>0& abs(D)>0.00005 & abs(D)<1.5 & abs(imag(D))<1e-6);
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.7 & abs(D)<1.3 & imag(D)>-1e-6 ); % find complex eigenvalues
end
if length(D)<=9
    h=1:length(D);
end
%h=1:min(9,length(D));
%% Draw Eigenfunctions
figure_num=9;
attachments=[];
for i=1:min(figure_num,length(h))
    %figure(i)
    %set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
    subplot(3,3,i)
    %tightsub(3,3,i,0.8)
    [str_nm,str_lamb,str_log,str_Terr,T]=titlestr(D(h(i)),a,b,m,n,md);
    for j=1
        %subplot(2,2,j)
        X=real(reshape(x_k,n,n));
        Y=real(reshape(y_k,n,n));
        if setup.leftU
            A=reshape(K*F(:,h(i)),n,n);
            %A=K*F(:,h(i));
        else
            A=reshape(F(:,h(i)),n,n);
            %A=F(:,h(i));
        end
        if j==1
            Z=real(A);
        elseif j==2
            Z=imag(A);
        elseif j==3
            Z=abs(A);
        elseif j==4
            Z=angle(A)/pi*180;
        end
        surf(X,Y,Z);
        %scatter3(x_k(:),y_k(:),Z(:),3,Z(:));%scatter3(X,Y,Z,S,C), 前三个参数是坐标，S是散点的size,C是颜色参数
        %hold on
        %Henon_attractors_draw([0 1 0],rate(Z,1),true);
%         Henon_T_draw(7,3,rate(Z,1.01));
%         Henon_T_linear_draw(7,3,rate(Z,1.02));
%         Henon_stable_manifold(7,3,rate(Z,1.03));
%         rz=inter(X,Y,Z,10,1);
        %Henon_iteration_forward_reverse;
        if j==1
            %ylabel('real')
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
%         xlim([-1.5,1.5])
%         ylim([-1.5,1.5])
        %axis([-1.5 1.5 -1.5 1.5])
        axis equal
        view(0,90)
        %view(-15,60)
        
    end
    %str='Eigenfunctions of Henon Map with Gauss Basis';
    % suptitle([str,' (n=',num2str(n),'^2,m=',num2str(m),'*2+1)']);
    %suptitle({[str,' (n=',num2str(n),',m=',num2str(m),'^2,d_j=',num2str(3/md,'%.4f'),')'];str_lamb});
    title(str_lamb);
    %xlabel(['T=',num2str(T)]);
    %saveas(gcf,['./temp/Henon_eigen_Gauss_n100m50md45_figure',num2str(i),'.png'])
end
% str='Eigenfunctions of Henon Map on Attractors with Gauss Basis';
str='Eigenfunctions and Stable Manifold of Henon Map';
% suptitle([str,' (n=',num2str(n),'^2,m=',num2str(m),'*2+1)']);
suptitle([str,' (n=',num2str(n),'^2,m=',num2str(m),'^2,d_j=',num2str(3/md,'%.4f'),') (T=7)']);
% saveas(gcf,['./temp/Henon_eigen_Gauss_manifold_n100m',num2str(m),'T7_3.png'])

function [str_nm,str_lamb,str_log,str_Terr,T]=titlestr(d,~,b,m,n,md)
d_abs=abs(d);
d_angle=angle(d)/pi*180;
N=log(b)/log(d_abs);
A=d_angle/N;
T=round(360/A);
err=abs((A*T-360))/360*100;
str_nm=['n=',num2str(n),'^2; m=',num2str(m),'^2; d_j=',num2str(3/md,'%.4f')];
str_lamb=['\lambda=',num2str(d_abs),'∠',num2str(d_angle),'°'];
str_log=['log_{',num2str(b),'}(',num2str(d_abs),')=',num2str(N),'; ',num2str(d_angle),'°/',num2str(N),'=',num2str(A),'°'];
str_Terr=['T=' num2str(360/A) '≈' num2str(T) '; err=' num2str(err) '%'];
end