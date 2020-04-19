function [K,L]=Koopman_basis_all(fun,param,opt)
x_k=param.x0(:);
x_l=fun(x_k);
n=param.n;
M=param.m;

count=1;
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
    if length(D)<=9
        h=1:length(D);
    end
    
    h=h(min(opt.all.choose,end));
%     if count==1
%         h=h(1);
%     else
%         h=h(min(opt.all.choose,end));
%     end
    
    param.h=h;
    param.D=D;
    param.figure_num=count;count=count+1;
    
    param.A=K*F(:,h);
    plotAll(param,opt);
end
suptitle(opt.title);

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


switch opt.all.deal
    case {'real'}
        A=real(A);
    case {'abs'}
        A=abs(A);
    case {'imag'}
        A=imag(A);
    case {'angle'}
        A=angle(A)*180/pi;
    otherwise
        error('options.all.deal param error.');
end

subplot(3,3,figure_num)
plot(x0,A);

d_abs=abs(D(h));
d_angle=angle(D(h))/pi*180;
str1=['m=',num2str(m(figure_num)),'; \lambda='];
str2=[num2str(d_abs),'б╧' num2str(d_angle),'бу'];
title([str1,str2]);

end

function saveit(param,opt)
h=gcf;
if opt.save.enabled
    temp=['_n',num2str(param.n),'_',toString(param.m,'-')];
    str=[opt.save.path,'/',opt.save.pre,temp,opt.save.suffix];
    %disp(str);
    saveas(h,str);
end
end

function str=toString(m,split)
str='m';
for i=1:length(m)-1
    str=[str,num2str(m(i)),split];
end
str=[str,num2str(m(end))];
end
