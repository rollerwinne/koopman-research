function [K,L]=Koopman_1d(fun,param,opt)
if length(param.m)>1
    param.multim.enabled=true;
end
if param.natural.enabled
    param.basis='natural';
    [K,L]=Koopman_natural(fun,param,opt);
    return
end
x_k=param.x0(:);
x_l=fun(x_k);
n=param.n;
M=param.m;
param.M=M;

cal=1;
figure;
set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
for m=M
    K=zeros(n,m);L=zeros(n,m);
    for j=1:m
        g=fun_basis(param.basis,j,m);
        K(:,j)=g(x_k);
        L(:,j)=g(x_l);
    end
    U=pinv(K)*L;
    [F,D]=eig(U);
    D=diag(D);
    h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6);
    maxh=opt.subp(1)*opt.subp(2);
    h=1:min(maxh,length(D));
    [~,idx]=sort(abs(D(h)),'descend');
    h=h(idx);
    
    if opt.multim.enabled
        H=h(min(opt.multim.choose,end));
    else
        H=h;
    end
    
    param.D=D;
    for h=H
        param.h=h;param.m=m;
        param.figure_num=cal;cal=cal+1;
        param.A=K*F(:,h);
        plotAll(param,opt);
    end
end
%suptitle(opt.title);

saveit(param,opt);
end

function [K,L]=Koopman_natural(fun,param,opt)
x0=param.natural.x0;
n=param.n;
M=param.m;param.M=M;
times=param.times;
dim=length(x0);

cal=1;
myfigure;
for m=M
    x_iter=param.natural.x0;
    K_x=zeros(dim*(m+n),1);L_x=K_x;K=zeros(dim*n,m);L=K;
    for i=1:m+n
        K_x((i-1)*dim+1:i*dim)=x_iter';
        x_temp=zeros(times,dim);
        for k=1:times %每个点演化迭代times次
            x_temp(k,:)=fun(x_iter);
        end
        x_iter=sum(x_temp,1)/times;%取平均
        L_x((i-1)*dim+1:i*dim)=x_iter';
    end
    for i=1:m
        K(:,i)=K_x( dim*(i-1)+1 : dim*(i-1)+1 + dim*n-1 );
        L(:,i)=L_x( dim*(i-1)+1 : dim*(i-1)+1 + dim*n-1 );
    end
    x_initial=K(:,1);param.x0=x_initial;
    
    U=pinv(K)*L;
    [F,D]=eig(U);
    D=diag(D);
    h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6);
    maxh=opt.subp(1)*opt.subp(2);
    h=1:min(maxh,length(D));
    [~,idx]=sort(abs(D(h)),'descend');
    h=h(idx);

    if opt.multim.enabled
        H=h(min(opt.multim.choose,end));
    else
        H=h;
    end
    
    param.D=D;
    for h=H
        param.h=h;param.m=m;
        param.figure_num=cal;cal=cal+1;
        param.A=K*F(:,h);
        plotAll(param,opt);
    end
end
%suptitle(opt.title);
saveit(param,opt);
end

function g=fun_basis(str,i,m)
switch str
    case {'Rectangle'}
        g = Rectangle(i,m);
    case {'Gauss'}
        g = Gauss(i,m);
    case {'Fourier'}
        g = Fourier(i,m);
    case {'Legendre'}
        g = Legendre(i,m);
    otherwise
        error('function name error');
end
end

function g = Rectangle(i,m)
if(i==m)
    g=@(x)(x>=(i-1)/m & x<=i/m)*sqrt(m);
else
    g=@(x)(x>=(i-1)/m & x<i/m)*sqrt(m);
end
end

function g = Gauss(i,m)
dj=1/m/2;
xj=linspace(1/2/m,1-1/2/m,m);
c=1/sqrt(dj*sqrt(pi));
g=@(x,m)c*exp(-(x-xj(i)).^2./(2.*dj.^2));
end

function g = Fourier(i,~)
if i==1
    g=@(x)ones(1,length(x))*sqrt(2)/2;
elseif mod(i,2)==0
    g=@(x)sqrt(2)*cos(2*pi*floor(i/2).*x);
else
    g=@(x)sqrt(2)*sin(2*pi*floor(i/2).*x);
end
end

function g = Legendre(i,~)
getFirst=@(x)x(1,:);
boundary=@(x)((x<0).*(x+1)+(x>1).*(x-1)+(x>0 & x<1).*x);
g=@(x)getFirst(legendre(i,boundary(x),'norm'));
end

function h=plotAll(param,opt)
x0=param.x0;
A=param.A;
h=param.h;
D=param.D;
n=param.n;
m=param.m;
figure_num=param.figure_num;


switch opt.multim.deal
    case {'real'}
        A=real(A);
    case {'abs'}
        A=abs(A);
    case {'imag'}
        A=imag(A);
    case {'angle'}
        A=angle(A)*180/pi;
    otherwise
        error('options.multim.deal param error.');
end

subplot(opt.subp(1),opt.subp(2),figure_num)

if param.natural.enabled
    [~,x0,A]=plotsort(x0,A,'k');
else
    plot(x0,A,'k','LineWidth',2);
end
if opt.peak.enabled
   [xp,yp]=draw_peaks(A,x0);
   [xt,yt]=draw_trough(A,x0);
end
if opt.boundary.enabled
    %D1_boundary_draw(opt.boundary.fun,1:min((m+1),9),(max(A)+min(A))/2,opt.boundary.color);
    bound=boundary_draw(opt.boundary.fun,1,min(m+1,6),(max(A)+min(A))/2);
end
% disp('-----m=')
% m
% disp('bound')
% bound(3:end)'
% [~,hi]=sort(bound(3:end));
% disp('peaks')
% temp=sort([xp;xt]);
% E(hi)=temp;
% E=E(:);
% E
% disp('-----end')

%xlabel('x');he=ylabel('$\phi\left(x\right)$');set(he,'Interpreter','latex');
set(gca,'yticklabel',num2str(get(gca,'YTick')','%.2f'));
d_abs=abs(D(h));
d_angle=angle(D(h))/pi*180;
str1=['m=',num2str(m),', \lambda='];
str2=[num2str(d_abs),'∠' num2str(d_angle),'°'];
title([str1,str2]);
sciformat(20)
end

function [h,x0,A]=plotsort(x0,A,color)
[x0,index]=sort(x0);
A=A(index);
h=plot(x0,A,color,'LineWidth',2);
end

function saveit(param,opt)
h=gcf;
temp=['_n',num2str(param.n),'_',toString(param.M,'-')];
str=[opt.save.pre,temp];
if opt.save.enabled
    savesci(str,h);
else
    disp(str);
end
end

function str=toString(m,split)
str='m';
for i=1:length(m)-1
    str=[str,num2str(m(i)),split];
end
str=[str,num2str(m(end))];
end

function [xt,yt]=draw_trough(A,x0)
hold on
[~,locs] = findpeaks(-A);
xt=x0(locs);
yt=A(locs);
plot(xt,yt,'bs');
end

function [xp,yp]=draw_peaks(A,x0)
hold on
[~,locs] = findpeaks(A);
xp=x0(locs);
yp=A(locs);
plot(xp,yp,'bs');
end
