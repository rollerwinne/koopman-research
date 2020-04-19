function [K,L]=Koopman_basis(fun,param,opt)
x_k=param.x0(:);
x_l=fun(x_k);
n=param.n;
m=param.m;

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
if length(D)<=9
    h=1:length(D);
end

param.h=h;
param.D=D;
for i=1:min(length(h),9)
    param.figure_num=i;
    param.A=K*F(:,h(i));
    plotDim(param,opt);
end

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
g=@(x)getFirst(legendre(i,x,'norm'));
end


function h=plotDim(param,opt)
x0=param.x0;
A=param.A;
dim=param.dim;
h=param.h;
D=param.D;
n=param.n;
m=param.m;
times=param.times;
figure_num=param.figure_num;

figure;
set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);

A_real=real(A);
A_abs=abs(A);

switch dim
    case 1
        subplot(dim,2,1)
        plot(x0,A_real);
        ylabel('x-real');
        subplot(dim,2,2)
        plot(x0,A_abs);
        ylabel('x-abs');
    case 2
        x_str=['x','y'];
        for i=1:dim
            subplot(dim,2,i*2-1)
            scatter3(x0(1:dim:end),x0(2:dim:end),A_real(i:dim:end),3,A_real(i:dim:end));
            ylabel([x_str(i),'-real']);
            graph_control_2D(opt);
            subplot(dim,2,i*2)
            scatter3(x0(1:dim:end),x0(2:dim:end),A_abs(i:dim:end),3,A_abs(i:dim:end));
            ylabel([x_str(i),'-abs']);
            graph_control_2D(opt);
        end
    case 3
        x_str=['x','y','z'];
        for i=1:dim
            subplot(2,dim,i)
            scatter3(x0(1:dim:end),x0(2:dim:end),x0(3:dim:end),3,A_real(i:dim:end));
            ylabel([x_str(i),'-real']);
            graph_control_3D(opt);
            subplot(2,dim,i+dim)
            scatter3(x0(1:dim:end),x0(2:dim:end),x0(3:dim:end),3,A_abs(i:dim:end));
            ylabel([x_str(i),'-abs']);
            graph_control_3D(opt);
        end
    otherwise
        for i=1:dim
            subplot(2,dim,i)
            plot(x0(i:dim:end),A_real(i:dim:end));
            ylabel(['x_',num2str(i),'-real']);
            subplot(2,dim,i+dim)
            plot(x0(i:dim:end),A_real(i:dim:end));
            ylabel(['x_',num2str(i),'-abs']);
        end
end

d_abs=abs(D(h(figure_num)));
d_angle=angle(D(h(figure_num)))/pi*180;
str1=['n=',num2str(n),'; m=',num2str(m),'; times=',num2str(times)];
str2=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
suptitle({str1;str2});

saveit(param,opt);
end

function saveit(param,opt)
h=gcf;
if opt.save.enabled
    str=[opt.save.path,'/',opt.save.pre,param_str(param),'figure',num2str(param.figure_num),opt.save.suffix];
    disp(str);
    %saveas(h,str);
end
end

function str=param_str(param)
str=['_n',num2str(param.n),'m',num2str(param.m),'_'];
end

function graph_control_2D(opt)
view(opt.view)
%view(options.view)
shading interp
colorbar
colormap(opt.map)
end

function graph_control_3D(opt)
colorbar
colormap(opt.map)
end