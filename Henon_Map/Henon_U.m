function [F,D,U] = Henon_U(x_k,x_l,y_k,y_l,m,md)
%% 求解U,及其特征值D和特征函数F
%   输入x_k：x_p
%   输入x_l: x_p+1

tic
% global md;
dj=3/md;    
% global m;    %基函数的个数

[xj,yj]=meshgrid(linspace(-1.5,1.5,m));%得到两个m*m的矩阵
mxy=length(x_k);
xj=xj(:);
yj=yj(:);

%% 高斯基
% figure(100)
for j=1:m*m
    for i=1:mxy
        K(i,j)=exp(-(x_k(i)-xj(j))^2/(dj^2)-(y_k(i)-yj(j))^2/(dj^2));        
    end
%     hold on
%     surf(reshape(x_k,sqrt(mxy),sqrt(mxy)),reshape(y_k,sqrt(mxy),sqrt(mxy)),reshape(K(:,j),sqrt(mxy),sqrt(mxy)))
%     pause(0.2)
%     shading interp
%     drawnow
end
% figure(101)
% surf(reshape(x_k,sqrt(mxy),sqrt(mxy)),reshape(y_k,sqrt(mxy),sqrt(mxy)),reshape(sum(K,2),sqrt(mxy),sqrt(mxy)))
% shading interp

for i=1:mxy
    for j=1:m*m
        L(i,j)=exp(-(x_l(i)-xj(j))^2/(dj^2)-(y_l(i)-yj(j))^2/(dj^2));        
    end
end
t=toc;
%disp(['Matric has already been caculated and cost ' num2str(t) ' seconds.'])

tic
iK=pinv(K);     %求K的广义逆：iK
U=L*iK;
[F,D]=eig(U);   %%求U的特征值D，特征函数F
D=diag(D);
%save(['.\data\Henon_Matrix_data_FD_n100m50md45a1.4b0.3.mat'],'F','D');
t=toc;
%disp(['Eigenfunctions has already been caculated and cost ' num2str(t) ' seconds.'])
end