clear;close all;clc

options=default_options;
[fun,param,options]=Logistic_map(options);

% param.n=1000; %演化格点维度
% param.m=3; %函数格点维度
% param.x0=[rand,rand]; %初始点
% param.times=1; %噪声平均次数
% options.view=[-15,60]; %画图角度
% options.map=jet; %色系


Koopman_basis_all(fun,param,options);

function [fun,param,options]=Tent_map(options)
alpha=1/2;
D=0;
fun=@(x)awgn((x<alpha)./alpha.*x+(x>=alpha).*(1/(1-alpha)-1/(1-alpha)*x),10*log10(1/D));

param.dim=1;
param.n=1000;
param.m=[2,3,4,8,10,16,20,50,100];
param.phase=[0,1];
param.x0=linspace(param.phase(1),param.phase(2),param.n);
param.times=1;
param.basis='Legendre';

options.title=['Eigenfunctions of Tent Map with ',param.basis,' Basis (n=',num2str(param.n),')'];
options.save.enabled=true;
options.save.path='./temp';
options.save.pre=['Tent_eigen_',param.basis];
options.save.suffix='.png';
end

function [fun,param,options]=Logistic_map(options)
alpha=4;
D=0;
fun=@(x)awgn(alpha.*x.*(1-x),10*log10(1/D));

param.dim=1;
param.n=1000;
param.m=[2,3,4,8,10,16,20,50,100];
param.phase=[0,1];
param.x0=linspace(param.phase(1),param.phase(2),param.n);
param.times=1;
param.basis='Legendre';

%options.title=['Eigenfunctions of Logistic Map with ',param.basis,' Basis (n=',num2str(param.n),')'];
options.title=['Eigenfunctions of Logistic Map with ',param.basis,' Basis (n=',num2str(param.n),',noise=',num2str(D),')'];
options.save.enabled=true;
options.save.path='./temp';
options.save.pre=['Logistic_eigen_noise_',param.basis];
options.save.suffix='.png';
end

function options=default_options
options.all.enabled=true;
options.all.deal='real';
options.all.choose=1;

options.save.enabled=false;
options.save.path='./temp';
options.save.suffix='.png';
end