function x=Lorenz_Poincare_next_npoints(x0,n)
% clear;clc
% x0=[-1,3,27];n=5;
rho=28;sigma=10;beta=8/3;
f=@(t,x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
tspan=[0,1000];
z0=x0(3);
options=odeset('Events',@(t,x)event_z(t,x,z0));
x_temp=x0;
x=[];
for i=1:n
	[~,~,~,Xend,~]=ode45(f,tspan,x_temp,options);
	x_temp=Xend(end,:);
	x=[x;x_temp];
end
end
function [value,isterminal,direction] = event_z(t,x,z0)
value = x(3)-z0;
isterminal=1;
direction =0;
end
