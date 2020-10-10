clear;close all;clc

options=default_options;
[fun,param,options]=Tent_map(options);
%Intermittency Tent_map Henon_map Henon3d_map
Koopman_natural(fun,param,options);

function [fun,param,options]=Intermittency(options)
fun=@(x)x./(1-x)*(x<=0.5&x>=0)+(1-x)/x*(x>0.5&x<=1);
param.n=1000;
param.m=4;
param.x0=rand;
param.times=1;
end

function [fun,param,options]=Tent_map(options)
alpha=1/2;
D=0;
fun=@(x)awgn((x<alpha)./alpha.*x+(x>=alpha).*(1/(1-alpha)-1/(1-alpha)*x),10*log10(1/D));

param.n=1000;
param.m=2;
param.x0=rand;
param.times=1;
end

function [fun,param,options]=Henon_map(options)
a=1.4;b=0.3;
D=0;
fun=@(x)[awgn(x(2)+1-a.*x(1).*x(1),10*log10(1/D));
    awgn(b*x(1),10*log10(1/D))];
param.n=5000;
param.m=2;
param.x0=[rand,rand];
param.times=1;

param.findpeak.enabled=true;
param.findpeak.nearby=25;
param.findpeak.color='black';

options.view=[-15,60];
options.map=jet;
options.save.enabled=false;
options.save.path='./temp';
options.save.pre='Henon_Koopman';
options.save.suffix='.png';
%view(-15,60)
%view(0,90)
end

function [fun,param,options]=Henon3d_map(options)
A=-1.86;
C=0.03;
D=0;
fun=@(x)[awgn(  x(2)  ,10*log10(1/D));
    awgn(  x(3)  ,10*log10(1/D));
    awgn(  0.72*x(1)+C*x(2)+A*x(3)-1.45*x(3).^2+0.515.*x(2).*x(3)-x(2).^2  ,10*log10(1/D))];
% A=0.82;
% C=2.06;
% fun2=@(x)[x(2);
%     x(3);
%     0.5.*x(1)+A*x(3)+C*x(2)-2.25*x(3).^3-2*x(2).^3];
% x0=[0.1,0,0];
param.n=5000;
param.m=2;
param.x0=[0.1,0,0];
param.times=1;
options.view=[0,0];
%options.view=[0,90];
options.map=jet;
options.save.enabled=true;
options.save.path='./temp';
options.save.pre='Henon3d2_Koopman';
options.save.suffix='.png';
%view(-15,60)
%view(0,90)
end

function options=default_options
param.findpeak.enable=false;
param.findpeak.nearby=100;
param.findpeak.color='black';

options.view=[-15,60];
options.map=jet;
options.save.enabled=false;
options.save.path='./temp';
options.save.pre='Henon3d_Koopman';
options.save.suffix='.png';
end