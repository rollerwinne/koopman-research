function [F,D]= Logistic_U_Gauss(x,m,md)
x_k=x(1,:);
x_l=x(end,:);
x0=linspace(0,1,m);

K=Logistic_G_Gauss(x_k,x0,m,md);
L=Logistic_G_Gauss(x_l,x0,m,md);

[F,D]=eig(L*pinv(K));
D=diag(D);
end

function G= Logistic_G_Gauss(x,xj,m,md)
dj=1/md;
G=zeros(length(x),m);
for i=1:length(x)
    for j=1:m
        G(i,j)=exp(-(x(i)-xj(j))^2/(2*dj^2));
    end
end
end


