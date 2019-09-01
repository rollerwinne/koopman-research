function [F,D,U,K,L]= Logistic_U_Gauss_sparse(x,m)
x_k=x(1,:);
x_l=x(end,:);
x0=linspace(1/2/m,1-1/2/m,m);

K=Logistic_G_Gauss(x_k,x0,m);
% figure(10)
% plot(x_k,K);
L=Logistic_G_Gauss(x_l,x0,m);
U=L*pinv(K);
[F,D]=eig(U);
D=diag(D);
end

function G= Logistic_G_Gauss(x,xj,m)
dj=1/m/2;
G=zeros(length(x),m);
for i=1:length(x)
    for j=1:m
        G(i,j)=exp(-(x(i)-xj(j))^2/(2*dj^2));
    end
end
end