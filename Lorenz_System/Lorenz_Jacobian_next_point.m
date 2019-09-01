function x=Lorenz_Jacobian_next_point(x0,z0,n)
% clear;clc
% x0=[-1,3,27];n=5;
if nargin<3;n=1;end
if nargin<2;z0=x0(3);end
tspan=[0,1000];
options=odeset('Events',@(t,x)event_z(t,x,z0));
x_temp=x0;
x=[];
for i=1:n
	[~,~,~,Xend,~]=ode45(@Lorenz_with_Jacobian,tspan,x_temp,options);
	x_temp=Xend(end,:);
	x=[x;x_temp];
end
end
function [value,isterminal,direction] = event_z(t,x,z0)
value = x(3)-z0;
isterminal=-1;
direction =1;
end