function [F,D,K]=Henon_basis_FDK(x,y,setup)
x_k=x(1,:);
y_k=y(1,:);%p时刻数据
x_l=x(end,:);
y_l=y(end,:);

switch setup.function
    case {'Poly'}
        % 需要最高次数
        [K,L] = Henon_KL_poly(x_k,x_l,y_k,y_l,setup.poly.power);
    case {'Gauss'}
        % 需要基函数个数及格点宽度指标
        [K,L] = Henon_KL_Gauss(x_k,x_l,y_k,y_l,setup.gauss.m,setup.gauss.md);
    case {'Fourier'}
        % 需要最高倍频数
        [K,L] = Henon_KL_Fourier(x_k,x_l,y_k,y_l,setup.fourier.m);
    case {'Legendre'}
        % 需要最高阶数
        [K,L]=Henon_KL_Legendre(x_k,x_l,y_k,y_l,setup.legendre.power);
    otherwise
        error('function not support');
end
if setup.leftU
    [F,D]=eig(pinv(K)*L); 
else
    [F,D]=eig(L*pinv(K));
end
D=diag(D);
end

%% Gauss Basis
function [K,L] = Henon_KL_Gauss(x_k,x_l,y_k,y_l,m,md)

dj=3/md;    
[xj,yj]=meshgrid(linspace(-1.5,1.5,m));%得到两个m*m的矩阵
xj=xj(:);
yj=yj(:);

for j=1:m*m
    K(:,j)=Gauss_basis_2d(x_k,y_k,xj(j),yj(j),dj);
    L(:,j)=Gauss_basis_2d(x_l,y_l,xj(j),yj(j),dj);
end
end

function y=Gauss_basis_2d(x,y,xj,yj,dj)
y=exp(...
( -(x-xj).^2-(y-yj).^2 )./(dj.^2)...
);
end

%% Fourier Basis
function [K,L] = Henon_KL_Fourier(x_k,x_l,y_k,y_l,m)

mxy=length(x_k);
[mx,ny]=meshgrid(1:m);
mx=mx(:);
ny=ny(:);
for j=1:m*m
    K_even(:,j)=Henon_cos(x_k,y_k,mx(j),ny(j));
    K_odd(:,j)=Henon_sin(x_k,y_k,mx(j),ny(j));
end
K=ones(mxy,2*m*m+1);
K(:,2:2:2*m*m)=K_even;
K(:,3:2:2*m*m+1)=K_odd;

for j=1:m*m
    L_even(:,j)=Henon_cos(x_l,y_l,mx(j),ny(j));
    L_odd(:,j)=Henon_sin(x_l,y_l,mx(j),ny(j));
end
L=ones(mxy,2*m*m+1);
L(:,2:2:2*m*m)=L_even;
L(:,3:2:2*m*m+1)=L_odd;
end

function res=Henon_cos(x,y,m,n)
res=cos(pi/1.5*(m*x+n*y));
end

function res=Henon_sin(x,y,m,n)
res=sin(pi/1.5*(m*x+n*y));
end

%% Ploynomial basis
function [K,L]=Henon_KL_poly(x_k,x_l,y_k,y_l,power)
K=Ploynomial_basis([x_k',y_k'],power);
L=Ploynomial_basis([x_l',y_l'],power);
end

function G=Ploynomial_basis(x,power)
%x=rand(4,1);dimention=3;power=2;
dimention=length(x(1,:));
n=length(x(:,1));
G=[];
for i=0:power
    if i==0
        G=[G,ones(n,1)];
    else
        P=Polynomial_fun(dimention,i,1);
        G_temp=[];
        for j=1:n
            G_temp=[G_temp;P(x(j,:))];
        end
        G=[G,G_temp];
    end
end
end

function [fun,l,f_split]=Polynomial_fun(d,p,flag) 
% d:dimention (d>=1)
% p;power (p>=1)
% flag=1(row) or flag=0(queue)
x_str=[];
x_addstr=[];
for i=1:d
    x_str=[x_str,' x',num2str(i)];
    if i==1
        x_addstr=[x_addstr,'x',num2str(i)];
    else
         x_addstr=[x_addstr,'+x',num2str(i)];
    end
end
str=['syms',x_str];
eval(str);
str=['f=expand((',x_addstr,')^',num2str(p),');'];
eval(str);
f_str=char(f);
f_str(isspace(f_str)) = [];
f_split=regexp(f_str,'+','split');
f_split=f_split';
f_funstr=[];
l=length(f_split);
for i=1:l
    if isstrprop(f_split{i}(1),'digit')
        f_where=strfind(f_split{i},'*');
        f_split{i}(1:f_where(1))=[];
    end
    f_split{i}=strrep(f_split{i},'*','.*');
    f_split{i}=strrep(f_split{i},'^','.^');
	[x_start,x_end]=regexp(f_split{i},'x[0-9]+');
    while ~isempty(x_start)
        f_split{i}=[f_split{i}(1:x_start(1)),'(',f_split{i}(x_start(1)+1:x_end(1)),')',f_split{i}(x_end(1)+1:end)];
        [x_start,x_end]=regexp(f_split{i},'x[0-9]+');
        %break;
    end
    if i==1
        f_funstr=f_split{i};
    elseif flag==1
        f_funstr=[f_funstr,',',f_split{i}]; % in a row
    else
        f_funstr=[f_funstr,';',f_split{i}]; % in a queue
    end
end
str=['fun=@(x)[',f_funstr,'];'];
eval(str);
end

%% Legendre basis
function [K,L]=Henon_KL_Legendre(x_k,x_l,y_k,y_l,power)
K=Legendre_basis_2d([x_k',y_k'],power);
L=Legendre_basis_2d([x_l',y_l'],power);
end

function G=Legendre_basis_2d(x,power)
width=[1.5,1.5];
center=[0,0];
n=length(x(:,1));
x_norm=zeros(n,2);
for i=1:2
    x_norm(:,i)=(x(:,i)-center(i))/width(i);
end
G=zeros(n,(power+2)*(power+1)/2);
count=1;
for i=0:power
    if i==0
        P=legendre(0,x_norm(:,1),'norm')/sqrt(width(1))/sqrt(width(2));
        %G=[G,P(1,:)'];
        G(:,count)=P(1,:)';count=count+1;
    else
        for nx=0:i
            PX=legendre(nx,x_norm(:,1),'norm')/sqrt(width(1));
            PY=legendre(i-nx,x_norm(:,2),'norm')/sqrt(width(2));
            G_temp=PX(1,:)'.*PY(1,:)';
            %G=[G,G_temp];
            G(:,count)=G_temp;count=count+1;
        end
    end
end
end