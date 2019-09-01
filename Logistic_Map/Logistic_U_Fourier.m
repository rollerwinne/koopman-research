function [F,D]= Logistic_U_Fourier(x,m)
x_k=x(1,:);
x_l=x(end,:);
K=Logistic_G_Fourier(x_k,m);
L=Logistic_G_Fourier(x_l,m);

[F,D]=eig(L*pinv(K));
D=diag(D);
end

function G= Logistic_G_Fourier(x,m)
G=zeros(length(x),m);
for j=1:m
    if j==1
        g=@(x)1;
    elseif mod(j,2)==0
        g=@(x)cos(floor(j/2).*x);
    else
        g=@(x)sin(floor(j/2).*x);
    end
    x=x(:);
    G(:,j)=g(x);
end
end