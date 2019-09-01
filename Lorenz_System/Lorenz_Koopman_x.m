function [T,X]=Lorenz_Koopman_x(tspan,x0)
rho=28;sigma=10;beta=8/3;
%f=@(x,y,z)[sigma*(y-x);x.*(rho-z)-y;x.*y-beta*z];
f=@(t,x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
%df=@(x,y,z)[-sigma,sigma,0;
%    rho-z,-1,-x;
%    y,x,-beta];
%x0=[-1,3,4];
[T,X]=rk_4(f,tspan,x0);
end