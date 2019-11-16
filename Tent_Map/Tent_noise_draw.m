%function [rou,X]=Tent_x_noise2(n,p,D)
% 测试数据
clear;close all
n=1000;%粒子个数
m=100;%基函数个数 6 100;100 10
times=10;%一次演化波包数
D=0.001;%噪声强度,信噪比的倒数


x0=linspace(0,1,n);
%[f,seq,sx]=Tents_function(7,D);           % Tents map
[f,seq,sx]=Tents_function_low(3,D);           % Tents map low
% f=@(x)awgn(1-2*abs(x-1/2),10*log10(1/D)); % Tent map with noise
% f=@(x)1-2*abs(x-1/2);                     % Tent map
% f=@(x)awgn(g(g(x)),10*log10(1/D));        % Tent map*2
% f=@(x)awgn(3.90.*x.*(1-x),10*log10(1/D)); % Logistic map with noise
% f=@(x)4.*x.*(1-x);                        % Logistic map
% f=@(x)awgn(2.5980762113533159402911695122588.*x.*(1-x).*(2-x),10*log10(1/D)); %偏移至0.41左右 with noise
% f=@(x)2.5980762113533159402911695122588.*x.*(1-x).*(2-x); %偏移至0.41左右
s=jet(n);
for m=[4,6,8,10,12,14,16,20,28,32,38,44,50,60,64,72,80,90,100]
    figure(m);
    subplot(3,3,1)
    plot(x0,f(x0))
    if exist('seq')
        hold on
        plot(seq,0,'r*')
    end
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
    U=U;
    [F,D]=eig(U);
    D=diag(D);
    h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6);
    if length(D)<=9
        h=1:length(D);
    end
    figure(m);
    for i=1:min(length(h),8)
        A=real(K*F(:,h(i)));
        subplot(3,3,i+1)
        % set(gcf,'outerposition',get(0,'screensize'));
        hh=plot(x0,A);
        hold on
        % draw_peaks(A,x0,seq,sx);
        % draw_trough(A,x0,seq,sx);
        [xp,yp]=draw_peaks(A,x0);
        [xt,yt]=draw_trough(A,x0);
        if (i==1)
            xpp=xp;
            ypp=yp;
            xtt=xt;
            ytt=yt;
            AA=A;
        end
        if exist('seq')
            plot(seq,min(A),'r*')
        end
        if exist('sx')
            draw_boundary(min(A),max(A),sx)
        end
        d_abs=abs(D(h(i)));
        d_angle=angle(D(h(i)))/pi*180;
        str1=['n=',num2str(n),'; m=',num2str(m)];
        str2=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
        title({str1;str2});
    end
    figure(101)
    subplot(121)
    hold on
    for i=1:length(ypp)
        y1index=floor( (ypp(i)-min(AA)) / (max(AA)-min(AA))*n+1);
        if y1index>n
            y1index=n;
        end
        plot(m,xpp(i),'*','Color',s(y1index,:))
    end
    subplot(122)
    hold on
    for i=1:length(ytt)
        y2index=floor( (ytt(i)-min(AA)) / (max(AA)-min(AA))*n+1);
        if y2index>n
            y2index=n;
        end
        plot(m,xtt(i),'*','Color',s(y2index,:))
    end
end
colormap(jet)
%figure(3);
%SVG_draw(U,K,n,m,0,0.1,1)

function [xt,yt]=draw_trough(A,x0,seq,sx)
[pks,locs] = findpeaks(-A);
xt=x0(locs);
yt=A(locs);
plot(xt,yt,'s');
hold on
if (nargin>2)
    [m,n]=size(sx);
    combine=[reshape(seq,length(seq),1);reshape(sx,m*n,1)];
    for i=1:length(xt)
        minindex=find(abs(abs(xt(i)-combine)) - min(abs(xt(i)-combine))<1e-10)
        diff=combine(minindex(1))-xt(i)
        text(xt(i),yt(i),num2str(abs(diff),'%5.4f'),'VerticalAlignment','bottom','rotation',90);
    end
end
end

function [xp,yp]=draw_peaks(A,x0,seq,sx)
[pks,locs] = findpeaks(A);
xp=x0(locs);
yp=A(locs);
plot(xp,yp,'s');
hold on
if (nargin>2)
    [m,n]=size(sx);
    combine=[reshape(seq,length(seq),1);reshape(sx,m*n,1)];
    for i=1:length(xp)
        minindex=find(abs(abs(xp(i)-combine)) - min(abs(xp(i)-combine))<1e-10)
        diff=combine(minindex(1))-xp(i)
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