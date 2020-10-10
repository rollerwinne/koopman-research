clear;close all
n=1000;%粒子个数
%m=50;%基函数个数 6 100;100 10
times=1;%一次演化波包数
d=0;%噪声强度,信噪比的倒数


x0=linspace(0,1,n);
%[f,seq,sx]=Tents_function(7,d);           % Tents map
%[f,seq,sx]=Tents_function(3,d);           % Tents map low()
f=@(x)awgn(1-2*abs(x-1/2),10*log10(1/d)); % Tent map with noise
% f=@(x)1-2*abs(x-1/2);                     % Tent map
%f=@(x)awgn(4.*x.*(1-x),10*log10(1/d));
S=load('tent_boundary_norepeat_x0.mat');
seq=[0,0.5,1];
sx=[S.X{4},S.X{4};S.X{5}];

figure(1)
plot(x0,f(x0));
if exist('seq')
    hold on
    plot(seq,0,'r*')
end

s=jet(n);YY=[];
M=[2,3,4,5,8,10,15,20];%[2,3,4,5,8,10,15,20]
for m=M
    figure(m);
    set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);

    X=zeros(n,times);
    K=zeros(n,m);
    L=zeros(n,m);
    for j=1:m %对于每个基函数
        % g_temp=Rect_fun(j,m);%rect基函数
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
    [fn1,fn2]=subfignum(min(length(h),9));
    figure(m);
    for i=1:min(length(h),9)
        A=real(K*F(:,h(i)));
        subplot(fn1,fn2,i)
        plot(x0,A,'k','LineWidth',2);
        hold on
        [xp,yp]=draw_peaks(A,x0);%画峰
        [xt,yt]=draw_trough(A,x0);%画谷
        %draw_boundary(min(A),max(A),sx);%画边界点
        plot(seq,min(A),'r*');
        %auto_level(f,[xp,xt],1e-2);
        auto_level(f,xp,yp,A,0.08,0.06*(max(A)-min(A)));
        auto_level(f,xt,yt,A,0.08,0.06*(max(A)-min(A)));
        
        d_abs=abs(D(h(i)));
        d_angle=angle(D(h(i)))/pi*180;
        str1=['m=',num2str(m),'; \lambda='];
        str2=[num2str(d_abs),'∠',num2str(d_angle),'°'];
        title([str1,str2]);
        sciformat(15)
    end
    saveas(gcf,['temp/Tent_auto_level_n1000_m',num2str(m),'.png'])
end
colormap(jet)

function hh=draw_pt(xpp,ypp,xtt,ytt,AA,s,n,m)
subplot(121)
hold on
for i=1:length(ypp)
    y1index=floor( (ypp(i)-min(AA)) / (max(AA)-min(AA))*n+1);
    if y1index>n
        y1index=n;
    end
    plot(m,xpp(i),'*','Color',s(y1index,:))
end
colorbar
title('波峰')
subplot(122)
hold on
for i=1:length(ytt)
    y2index=floor( (ytt(i)-min(AA)) / (max(AA)-min(AA))*n+1);
    if y2index>n
        y2index=n;
    end
    hh=plot(m,xtt(i),'*','Color',s(y2index,:));
end
colorbar
title('波谷')
end

function [xt,yt]=draw_trough(A,x0,fsize,seq,sx)
[pks,locs] = findpeaks(-A);
xt=x0(locs);
yt=A(locs);
plot(xt,yt,'s');
hold on
if (nargin>3)
    [m,n]=size(sx);
    combine=[reshape(seq,length(seq),1);reshape(sx,m*n,1)];
    for i=1:length(xt)
        minindex=find(abs(abs(xt(i)-combine)) - min(abs(xt(i)-combine))<1e-10)
        diff=combine(minindex(1))-xt(i)
        text(xt(i),yt(i),num2str(abs(diff),'%5.4f'),'FontSize',fsize,'VerticalAlignment','bottom','rotation',90);
    end
end
end

function [xp,yp]=draw_peaks(A,x0,fsize,seq,sx)
[pks,locs] = findpeaks(A);
xp=x0(locs);
yp=A(locs);
plot(xp,yp,'s');
hold on
if (nargin>3)
    [m,n]=size(sx);
    combine=[reshape(seq,length(seq),1);reshape(sx,m*n,1)];
    for i=1:length(xp)
        minindex=find(abs(abs(xp(i)-combine)) - min(abs(xp(i)-combine))<1e-10)
        diff=combine(minindex(1))-xp(i)
        text(xp(i),yp(i),num2str(abs(diff),'%5.4f'),'FontSize',fsize,'VerticalAlignment','bottom','rotation',90);
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
% x is a number
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