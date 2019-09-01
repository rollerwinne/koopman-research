function dF=Lorenz_with_Jacobian(t,x)
rho=28;sigma=10;beta=8/3;
F=@(x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
Df=@(x)[-sigma,sigma,0;
    rho-x(3),-1,-x(1);
    x(2),x(1),-beta];
DJ=@(J)[J(1),J(4),J(7);
    J(2),J(5),J(8);
    J(3),J(6),J(9)];
x_Lorenz=x(1:3);
x_J=x(end-8:end);
df_Lorenz=F(x_Lorenz);
df_Jacobian=Df(x_Lorenz)*DJ(x_J);
dF=[df_Lorenz(:);df_Jacobian(:)];
end

