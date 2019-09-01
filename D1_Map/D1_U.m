function [F,D] = D1_U(x_k,x_l,m,md)

dj=4/md;    
xj=linspace(-2,2,m);%得到两个m*m的矩阵
mx=length(x_k);
xj=xj(:);

for i=1:mx
    for j=1:m
        K(i,j)=exp(-abs(x_k(i)-xj(j))^2/(dj^2));        
    end
end

for i=1:mx
    for j=1:m
        L(i,j)=exp(-abs(x_l(i)-xj(j))^2/(dj^2));        
    end
end

[F,D]=eig(L*pinv(K));   %%求U的特征值D，特征函数F
D=diag(D);
end