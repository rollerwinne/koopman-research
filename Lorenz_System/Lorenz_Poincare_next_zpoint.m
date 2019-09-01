function x=Lorenz_Poincare_next_zpoint(x0,z0)
% clear;clc
% x0=[-1,3,27];
rho=28;sigma=10;beta=8/3;
f=@(t,x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
tspan=[0,100];
options=odeset('Events',@(t,x)event_z(t,x,z0));
[~,~,~,Xend,~]=ode45(f,tspan,x0,options);
% [T,X,Tend,Xend,evennum]=ode45(f,tspan,x0,options);
x=Xend(end,:);
end
function [value,isterminal,direction] = event_z(t,x,z0)
value = x(3)-z0;
isterminal=1;
direction =0;
end
