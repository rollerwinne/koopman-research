clear;close all
n=1000;%粒子个数
times=1;%一次演化波包数
D=0.000;%噪声强度,信噪比的倒数

x0=linspace(0,1,n);
%[f,seq,sx]=Tents_function(7,D);           % Tents map
[f,seq,sx]=Tents_function(3,D);           % Tents map low()
% f=@(x)awgn(1-2*abs(x-1/2),10*log10(1/D)); % Tent map with noise
% f=@(x)1-2*abs(x-1/2);                     % Tent map
% f=@(x)awgn(g(g(x)),10*log10(1/D));        % Tent map*2
% f=@(x)awgn(3.90.*x.*(1-x),10*log10(1/D)); % Logistic map with noise
% f=@(x)4.*x.*(1-x);                        % Logistic map
% f=@(x)awgn(2.5980762113533159402911695122588.*x.*(1-x).*(2-x),10*log10(1/D)); %偏移至0.41左右 with noise
% f=@(x)2.5980762113533159402911695122588.*x.*(1-x).*(2-x); %偏移至0.41左右
s=jet(n);
flag=0;
for m=[64,128]%[4,6,8,10,12,14,16,20,28,32,38,44,50,60,64,72,80,90,100]
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
    h=find(abs(D)>0.001 & abs(D)<1.3 & abs(imag(D))<1e-2);
    if length(D)<9
        h=1:length(D);
    end
    if flag==0
        flag=1;
        m_pre=m;
        for j=1:min(length(h),9)
            choose=j;
            X_eigen(:,choose)=real(K*F(:,h(choose)));
            X_para(:,choose)=real(F(:,h(choose)));
            lambda(choose)=D(h(choose));
        end
    else
        for j=1:min(length(lambda),9)
            subplot(3,3,j)
            plot(x0,X_eigen(:,j))
            hold on
            hh=Koopman_findeigen_draw(K,U,X_eigen(:,j),X_para(:,j),n,lambda(j),0.01,0);
            %legend(['by m=',num2str(m_pre)],['finding m=',num2str(m)])
        end
        suptitle(['finding m=',num2str(m),' by m=',num2str(m_pre)]);
        set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
        saveas(hh,['temp4/Tent_findeigen_m',num2str(m_pre),'m',num2str(m),'.png']);
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