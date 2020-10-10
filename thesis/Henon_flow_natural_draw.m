%% Initialization
clc;close all;clear
%clc;clear
tic;timestart=char(datetime('now'));
%disp('The running program is from ZC. 么么哒')
%% Parameter settings
a=1.4;b=0.3;n=10000;m=3;itern=11;
times=1;d=0;
Iter=m-1;
s=hsv(Iter);

Attr=load('./data/Henon_attractors_data_xy.mat'); % 吸引子数据载入
x_attr=Attr.x;y_attr=Attr.y;
Peri=load('./data/Henon_period_orbrits_P_1_0.3_-1_-0.3.mat');
P=Peri.P;
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;

K=zeros(2*n,m);L=K;
for j=1:m
    K(:,j)=Henon_attractors_flow(n,itern+j-1);
    L(:,j)=Henon_attractors_flow(n,itern+j);
end



%% Caculate Eigenfunction
%[F,D,U] = Henon_U(x_k,x_l,y_k,y_l,m,md);
%save('.\data\Henon_Matrix_data_on_attractors_FD_n50m50.mat','F','D');
%load('.\data\Henon_Matrix_data_on_attractors_FD_n50m50.mat');
U=pinv(K)*L;
[F,D]=eig(U);
D=diag(D);
%h=1:length(D);
%h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6);


%% Data processing
choose='complex';
if strcmp(choose,'real')==1
    %h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
    h=find(real(D)>0& abs(D)>0.0001 & abs(D)<1.2 & abs(imag(D))<1e-6);
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.5 & abs(D)<1.05 & imag(D)>-1e-6 ); % find complex eigenvalues
end
if length(D)<=9
    h=1:length(D);
end
[~,idx]=sort(abs(D(h)),'descend');
h=h(idx);
%% Draw Eigenfunctions
[fn1,fn2]=subfignum(length(h));
figure_num=9;
%tightsub(3,3,m,0.8);
%set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
for i=1:min(figure_num,length(h))
    %subplot(fn1,fn2,i)
    myfigure;
    for j=2
        X=K(1:2:end,1);
        Y=K(2:2:end,1);
        Z_temp=K*F(:,h(i));
        if j==1
            Z=real(Z_temp(1:2:end,1));
        elseif j==2
            Z=real(Z_temp(2:2:end,1));
        elseif j==3
            Z=real(Z_temp(1:2:end,1));
        elseif j==4
            Z=real(Z_temp(2:2:end,1));
        end
        scatter3(X(:),Y(:),Z(:),3,Z(:));
        z_min=min(Z(:));z_max=max(Z(:));
        R=Henon_flowpeaks_draw(X(:),Y(:),Z(:),'black');
        hold on
        plot3([Attr.x1;Attr.x2],[Attr.y1;Attr.y2],rate(Z,0.01)*ones(2,1),'r*');

        for iter=1:Iter
            disp([num2str(i),':',num2str(-iter+1)])
            [Q,~]=Henon_boundary_draw(-iter+1,rate(Z,0.01),rate(Z,1.05),s(iter,:));
            res=find_nearest(Q,R);
            disp(res);
        end
        

        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        if j==1
            scatter3(x_attr(:),y_attr(:),z_min*ones(length(x_attr(:)),1),3,min(Z(:))*ones(length(x_attr(:)),1));%画吸引子
            view(-15,60)
        elseif j==2
            scatter3(x_attr(:),y_attr(:),z_min*ones(length(x_attr(:)),1),3,min(Z(:))*ones(length(x_attr(:)),1));%画吸引子
            view(-15,60)
        elseif j==3
            view(0,90)
        elseif j==4
            view(0,90)
        end
        d_abs=abs(D(h(i)));
        d_angle=angle(D(h(i)))/pi*180;
        title(['m=',num2str(m),' \lambda=',num2str(d_abs),'∠' num2str(d_angle),'°']);
        shading interp
        colorbar
        colormap(jet)
    end

end
str=['Eigenfunction of Henon Map with Natural Basis (n=',num2str(n),',m=',num2str(m),')'];
%str=['Eigenfunction and Boundarys of Henon Map with Natural Basis (n=',num2str(n),')'];
%suptitle(str)
str=['./temp/Henon_eigen_natural_attr_n',num2str(n),'m',num2str(m)];
%str=['./temp/Henon_eigen_natural_attr_n',num2str(n)];
%saveas(gcf,[str,'.png']);

function maketitle(d,~,b,m,n)
d_abs=abs(d);
d_angle=angle(d)/pi*180;
N=log(b)/log(d_abs);
A=d_angle/N;
T=round(360/A);
err=abs((A*T-360))/360*100;
str0=['n=',num2str(n),'; m=',num2str(m)];
str1=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
str2=['log' '_{' num2str(b) '}(' num2str(d_abs) ')=' num2str(N) '; ' num2str(d_angle) '°/' num2str(N) '=' num2str(A) '°'];
str3=['T=' num2str(360/A) '≈' num2str(T) '; err=' num2str(err) '%'];
%suptitle({str0;[str2,'; ',str3];str1})
str=['Eigenfunction of Henon Map with Natural Basis (n=',num2str(n),',m=',num2str(m),')'];
suptitle({str;str1})
end

function z=rate(Z,r)
zmin=min(Z(:));
zmax=max(Z(:));
z=zmin+r*(zmax-zmin);
end