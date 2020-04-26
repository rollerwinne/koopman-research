function [F,D,K,x_k]=Lorenz_natural_U(x,n,m)
x0=[-1,3,4];
dim=length(x0);

K_x=zeros(dim*(m+n),1);L_x=K_x;K=zeros(dim*n,m);L=K;
for i=1:m+n
    K_x((i-1)*dim+1:i*dim)=x(i,:)';
    L_x((i-1)*dim+1:i*dim)=x(i+1,:)';
end
for i=1:m
    K(:,i)=K_x( dim*(i-1)+1 : dim*(i-1)+1 + dim*n-1 );
    L(:,i)=L_x( dim*(i-1)+1 : dim*(i-1)+1 + dim*n-1 );
end
x_k=K(:,1);

U=pinv(K)*L;
[F,D]=eig(U);
D=diag(D);
end