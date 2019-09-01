function [F,D,x_function_length]=Lorenz_Koopman_U_legendre(x_k,x_l,power)
K=Legendre_basis_3d(x_k,power);
L=Legendre_basis_3d(x_l,power);
x_function_length=length(K(1,:));
%% Eigenvalues and Eigenfuncitons of Koopman Operator
iK=pinv(K);     %求K的广义逆：iK
U=L*iK;
[F,D]=eig(U);   %%求U的特征值D，特征函数F
D=diag(D);
end