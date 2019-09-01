clear;close all
n=1000;p=1;m=100;
xx=Tent_x(n,p);
% xx=Logistic_x(4,n,p);
%[~,~,U,K,~,x_marker] = Tent_U_Rectangle_leftU(xx,m);
[~,~,U,K,~,x_marker] = Tent_U_Gauss_leftU(xx,m);

A=U';
[F,D]=eig(A);
D=diag(D);
h=find( abs(D)<1.3 & imag(D)>=0);
x0=linspace(0,1,n);
figure(1)
for i=1:min(m,16)
    subplot(4,4,i)
    pA=real(K*F(:,h(i)));
    plot(x0,pA)
    hold on
    plot([0,1],[0,0],'r')
    title(['\lambda=',num2str(abs(D(i))),'∠',num2str(angle(D(i))/pi*180),'°'])
end

% figure
% spy(A)
% title('U稀疏矩阵的非零值分布')