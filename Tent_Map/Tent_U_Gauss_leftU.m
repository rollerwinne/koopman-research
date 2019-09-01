function [F,D,U,K,L,x0]= Tent_U_Gauss_leftU(x,m)
x_k=x(1,:);
x_l=x(end,:);
x0=linspace(1/2/m,1-1/2/m,m);

K=Tent_G_Gauss(x_k,x0,m);
figure(10)
% set(gcf,'outerposition',get(0,'screensize'));
hhh=plot(x_k,K);
title(['Gauss Basis Function of m=',num2str(m)])
% str=['.\temp\Logistic_Gauss_basis_m',num2str(m)];
% saveas(hhh,[str,'.fig'])
% saveas(hhh,[str,'.png'])
L=Tent_G_Gauss(x_l,x0,m);
U=pinv(K)*L;
[F,D]=eig(U);
D=diag(D);
end

function G= Tent_G_Gauss(x,xj,m)
dj=1/m/2;
G=zeros(length(x),m);
for i=1:length(x)
    for j=1:m
        G(i,j)=exp(-(x(i)-xj(j))^2/(2*dj^2));
    end
end
end