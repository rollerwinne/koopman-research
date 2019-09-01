function [F,D,x_function_length,dj]=Lorenz_Koopman_U_natural(x_k,x_l)
tic
rho=28;sigma=10;beta=8/3;
%f=@(x,y,z)[sigma*(y-x);x.*(rho-z)-y;x.*y-beta*z];
f=@(t,x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
x0=[-1,3,4];
t=30;
[~,X]=ode45(f,[0,t],x0);
X_diff=diff(X);
dj=mean(sqrt(X_diff(:,1).^2+X_diff(:,2).^2+X_diff(:,3).^2));
for i=1:3
    x_function_lattice{i}=X(:,i);
end
x_function_length=length(x_function_lattice{1});
x_evolution_length=length(x_k);

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