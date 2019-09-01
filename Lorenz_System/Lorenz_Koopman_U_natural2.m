function [F,D,x_function_length,dj]=Lorenz_Koopman_U_natural2(x_k,x_l)
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
dj=sqrt(X_diff(:,1).^2+X_diff(:,2).^2+X_diff(:,3).^2);
x_function_lattice=X;

x_function_length=size(x_function_lattice,1)-1;
x_evolution_length=length(x_k);

%% Gauss Basic function
K=[];L=[];
for j=1:x_function_length
    K=[K,exp( (-(x_k(:,1)-x_function_lattice(j,1)).^2 ...
        -(x_k(:,2)-x_function_lattice(j,2)).^2 ...
        -(x_k(:,3)-x_function_lattice(j,3)).^2 ) /(dj(j)^2) )];
    L=[L,exp( (-(x_l(:,1)-x_function_lattice(j,1)).^2 ...
        -(x_l(:,2)-x_function_lattice(j,2)).^2 ...
        -(x_l(:,3)-x_function_lattice(j,3)).^2 ) /(dj(j)^2) )];
end

% for i=1:x_evolution_length
%     for j=1:x_function_length
%         K(i,j)=exp( (-(x_k(i,1)-x_function_lattice(j,1))^2 ...
%             -(x_k(i,2)-x_function_lattice(j,2))^2 ...
%             -(x_k(i,3)-x_function_lattice(j,3))^2 ) /(dj(j)^2) );
%         L(i,j)=exp( (-(x_l(i,1)-x_function_lattice(j,1))^2 ...
%             -(x_l(i,2)-x_function_lattice(j,2))^2 ...
%             -(x_l(i,3)-x_function_lattice(j,3))^2 ) /(dj(j)^2) );
%     end
% end

t=toc;
disp(['Lorenz_Koopman_U: Matric has been caculated and cost ' num2str(t) ' seconds.'])
%% Eigenvalues and Eigenfuncitons of Koopman Operator
%U=L/K;
%U=L*pinv(K);
[F,D]=eig(L/K);
D=diag(D);
t=toc;
disp(['Lorenz_Koopman_U: Eigenfunctions has been caculated and cost ' num2str(t) ' seconds.'])
end