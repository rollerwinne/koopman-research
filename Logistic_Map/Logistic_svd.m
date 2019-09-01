clear;clc;
n=1000;m=100;alpha=4;
x0=rand(1,1);  
f=@(x)alpha.*x.*(1-x);

x=x0;
X=x0;
for i=1:m+n-1
    x=f(x);
    X=[X;x];
end

M=[];
for i=1:m
    M(:,i)=X(i:i+n-1);
end

[U,S,V]=svd(M);
norm(M-U*S*V')%比较奇异值分解误差
select=10;
K=V(:,1:select)';
L=V(:,2:select+1)';
U=L*pinv(K);
[F,D]=eig(U);

% Logistic_U_Gauss_sparse(x,m)
% V1=combine(V);
% V2=split(V1,m);
% [F,D,U,K,L]= Logistic_U_Gauss_sparse(V1',10);
% temp=abs(F(:,1));

function new_V=combine(V)
[n,m]=size(V);
new_V(:,1)=reshape(V(:,1:2:end),n*m/2,1);
new_V(:,2)=reshape(V(:,2:2:end),n*m/2,1);
end

function new_V=split(V,n)
[l,m]=size(V);
s=l/n;
new_V=zeros(n,2*s);
new_V(:,1:2:end)=reshape(V(:,1),n,s);
new_V(:,2:2:end)=reshape(V(:,2),n,s);
end