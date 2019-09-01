function [F,D,U,x_function_num]=Lorenz_Koopman_U(x_k,x_l,dj)
tic
x_boundary=[-20,20;-30,30;0,50];
x_start=x_boundary(:,1);
x_function_num=ceil((x_boundary(:,2)-x_boundary(:,1))./dj);
x_function_num_temp=[];
for i=1:3
    x_function_num_temp{i}=x_start(i):dj:x_start(i)+dj*(x_function_num(i)-1);
end
[x_function_lattice{1},x_function_lattice{2},x_function_lattice{3}]= ...
    meshgrid(x_function_num_temp{1},x_function_num_temp{2},x_function_num_temp{3});
for i=1:3
    x_function_lattice{i}=x_function_lattice{i}(:);
end
x_function_length=length(x_function_lattice{1});
x_evolution_length=length(x_k);
% x_boundary_aver=mean(x_boundary(:,2)-x_boundary(:,1));
% dj=x_boundary_aver/md;
% x_function_num=ceil((x_boundary(:,2)-x_boundary(:,1))/dj);
% [x_function_lattice{1},x_function_lattice{2},x_function_lattice{3}]= ...
%     meshgrid(linspace(x_boundary(1,1),x_boundary(1,2),x_function_num(1)), ...
%     linspace(x_boundary(2,1),x_boundary(2,2),x_function_num(2)), ...
%     linspace(x_boundary(3,1),x_boundary(3,2),x_function_num(3)) );
% for i=1:3
%     x_function_lattice{i}=x_function_lattice{i}(:);
% end
% x_function_length=length(x_function_lattice{1});
% x_evolution_length=length(x_k);

%% Gauss Basic function
K=[];L=[];
for i=1:x_evolution_length
    for j=1:x_function_length
        K(i,j)=exp( (-(x_k(i,1)-x_function_lattice{1}(j))^2 ...
            -(x_k(i,2)-x_function_lattice{2}(j))^2 ...
            -(x_k(i,3)-x_function_lattice{3}(j))^2 ) /(dj^2) );
    end
end

for i=1:x_evolution_length
    for j=1:x_function_length
        L(i,j)=exp( (-(x_l(i,1)-x_function_lattice{1}(j))^2 ...
            -(x_l(i,2)-x_function_lattice{2}(j))^2 ...
            -(x_l(i,3)-x_function_lattice{3}(j))^2 ) /(dj^2) );        
    end
end
t=toc;
disp(['Lorenz_Koopman_U: Matric has been caculated and cost ' num2str(t) ' seconds.'])
%% Eigenvalues and Eigenfuncitons of Koopman Operator
iK=pinv(K);     %求K的广义逆：iK
U=L*iK;
[F,D]=eig(U);   %%求U的特征值D，特征函数F
D=diag(D);
t=toc;
disp(['Lorenz_Koopman_U: Eigenfunctions has been caculated and cost ' num2str(t) ' seconds.'])
end