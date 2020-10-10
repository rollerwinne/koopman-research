clear;close all;clc%噪声情况下，映射在不同基的本征函数，及峰谷位置
n=1000;%粒子个数
%m=100;%基函数个数 6 100;100 10
times=1;%一次演化波包数
d=0;%噪声强度,信噪比的倒数

x0=linspace(0,1,n);
%[f,seq,sx]=Tents_function(5,d);           % Tents map
%[f,seq,sx]=Tents_function_low(5,d);           % Tents map low()
%f=@(x)awgn(1-2*abs(x-1/2),10*log10(1/d)); % Tent map with noise
% f=@(x)1-2*abs(x-1/2);                     % Tent map
%f=@(x)awgn(g(g(x)),10*log10(1/d));        % Tent map*2
f=@(x)awgn(4.*x.*(1-x),10*log10(1/d)); % Logistic map with noise
%f=@(x)4.*x.*(1-x);                       % Logistic map
S=load('logistic_boundary_x0.mat');
seq=S.X{5};
% S2=load('logistic_boundary_x0.75.mat');
% seq2=S2.X{4};

s=jet(n);YY=[];
M=[2,4,8,16];%,14,16,20,28,32,38,44,50,60,64,72,80,90,100];
%M=[2,4,10,20];
figure;set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
plot(x0,f(x0))
if exist('seq')
    hold on
    plot(seq,0.5,'r*');
    %plot(seq2,1,'g*');
    %draw_boundary(0,1,sx)
    %title('Phase Space of Entire-Half Map');
end
%saveas(gcf,['./temp/Tents5_phase_d',num2str(d),'.png']);
cal=1;
figure;
for m=M
    subplot(2,2,cal);cal=cal+1;
    %figure(m);
    set(gcf,'outerposition',get(0,'screensize')-[0,0,1440-900,0]);
    X=zeros(n,times);
    K=zeros(n,m);
    L=zeros(n,m);
    for j=1:m %对于每个基函数
        %     g_temp=Rect_fun(j,m);%rect基函数
        g_temp=Gauss_fun(m);%Gauss基函数
        for i=1:length(x0) %对于每个相空间的点
            K(i,j)=g_temp(x0(i),j,m);
            x_temp=[];
            for k=1:times %每个点演化迭代times次
                x_temp=[x_temp,g_temp(f(x0(i)),j,m)];
            end
            %x_temp=x_temp(x_temp<=1 & x_temp>=0);
            L(i,j)=sum(x_temp)/times;
        end
    end
    U=pinv(K)*L;
    [F,D]=eig(U);
    D=diag(D);
    h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6);
    %h=find(abs(D)>0.001 & abs(D)<1.3 & abs(imag(D))<1e-13);
    h=1:min(9,length(D));
    [~,idx]=sort(D(h),'descend');
    h=h(idx);
    %[fn1,fn2]=subfignum(min(length(h),9));
    %figure(m);
    for i=1%1:min(length(h),9)
        A=real(K*F(:,h(i)));
        %subplot(fn1,fn2,i)
        hh=plot(x0,A,'k');
        hold on
        [xp,yp]=draw_peaks(A,x0);
        [xt,yt]=draw_trough(A,x0);
        
        if exist('seq')%边界点
%           plot(seq,(min(A)+max(A))/2,'r*')
            boundary_draw('logistic',1,5,(min(A)+max(A))/2);
            %plot(seq,min(A),'r*')
        end
        if exist('seq2')%边界点
            plot(seq2,max(A),'g*')
        end
        if exist('sx')%边界点的原像点
            draw_boundary(min(A),max(A),sx)
        end
        xlabel('x');he=ylabel('$\phi\left(x\right)$');set(he,'Interpreter','latex');
        d_abs=abs(D(h(i)));
        d_angle=angle(D(h(i)))/pi*180;
        str1=['m=',num2str(m)];
        str2=[', \lambda=',num2str(d_abs,'%.4f'),'∠' num2str(d_angle),'°'];
        title([str1,str2]);
    end
    %suptitle(['Eigenfunctions of Tents Map with noise (n=',num2str(n),',noise=',num2str(d),')']);
    %suptitle(['Eigenfunctions of Tents Map (n=',num2str(n),')']);
    %filename=['Tents5_eigen_noise_n1000m',num2str(m),'d',num2str(d),'.png'];
    %saveas(hh,['./temp/',filename]);
    sciformat(20);
end
%colormap(jet)
%suptitle(['Eigenfunctions of Tents Map with noise (n=',num2str(n),',noise=',num2str(d),')']);
%suptitle(['Eigenfunctions of Tents Map (n=',num2str(n),')']);
% filename=['Tents5_eigen_noise_n1000_m',num2str(m),'d',num2str(d),'.png'];
filename=['Logistic_eigen_Gauss_n1000_m',toString(M,'-'),'_d',num2str(d)];
%savesci(filename);

function [xt,yt]=draw_trough(A,x0,seq,sx)
[~,locs] = findpeaks(-A);
xt=x0(locs);
yt=A(locs);
plot(xt,yt,'bs');
hold on
if (nargin>2)
    [m,n]=size(sx);
    combine=[reshape(seq,length(seq),1);reshape(sx,m*n,1)];
    for i=1:length(xt)
        minindex=find(abs(abs(xt(i)-combine)) - min(abs(xt(i)-combine))<1e-10);
        diff=combine(minindex(1))-xt(i);
        text(xt(i),yt(i),num2str(abs(diff),'%5.4f'),'VerticalAlignment','bottom','rotation',90);
    end
end
end

function [xp,yp]=draw_peaks(A,x0,seq,sx)
[~,locs] = findpeaks(A);
xp=x0(locs);
yp=A(locs);
plot(xp,yp,'bs');
hold on
if (nargin>2)
    [m,n]=size(sx);
    combine=[reshape(seq,length(seq),1);reshape(sx,m*n,1)];
    for i=1:length(xp)
        minindex=find(abs(abs(xp(i)-combine)) - min(abs(xp(i)-combine))<1e-10);
        diff=combine(minindex(1))-xp(i);
        text(xp(i),yp(i),num2str(abs(diff),'%5.4f'),'VerticalAlignment','bottom','rotation',90);
    end
end
end

function draw_boundary(min,max,sx)
l=length(sx(:,1));
s=cool(l);
for i=1:l
    plot(sx(i,:),(max-min)/l*i+min,'*','color',s(i,:));
end
end

function g = Rect_fun(i,m)
if(i==m)
    g=@(x,i,m)(x>=(i-1)/m & x<=i/m)*sqrt(m); %x坐标，第i个基函数
else
    g=@(x,i,m)(x>=(i-1)/m & x<i/m)*sqrt(m);
end
end

function g= Gauss_fun(m)
dj=1/m/2;
xj=linspace(1/2/m,1-1/2/m,m);
g=@(x,i,m)exp(-(x-xj(i)).^2./(2.*dj.^2));
end