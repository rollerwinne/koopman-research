function [F,D,x_function_length]=Lorenz_Koopman_U_poly(x_k,x_l,power)
K=Ploynomial_basis(x_k,power);
L=Ploynomial_basis(x_l,power);
x_function_length=length(K(1,:));
%% Eigenvalues and Eigenfuncitons of Koopman Operator
iK=pinv(K);     %求K的广义逆：iK
U=L*iK;
[F,D]=eig(U);   %%求U的特征值D，特征函数F
D=diag(D);
end