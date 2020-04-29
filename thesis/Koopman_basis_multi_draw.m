clear;close all;clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=default_options;                    %
[fun,param,options]=Logistic_map(options); %
Koopman_basis_multi(fun,param,options);     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tent map
function [fun,param,options]=Tent_map(options)
alpha=1/2;
D=0;
fun=@(x)awgn((x<alpha)./alpha.*x+(x>=alpha).*(1/(1-alpha)-1/(1-alpha)*x),10*log10(1/D));

param.dim=1;
param.n=1000;
param.m=[1,2,3,4,5,6,7,8,9];%[2,3,4,5,8,10,15,20,50];%[2,3,4,8,10,16,20,50,100];
param.phase=[0,1];
param.x0=linspace(param.phase(1),param.phase(2),param.n);
param.times=1;
param.basis='natural';%'natural'

param.natural.enabled=true;%自然基
param.natural.x0=rand;

options.multim.enabled=true;%支持不同的基函数数量
options.multim.deal='real';
options.multim.choose=1;

options.boundary.enabled=true;
options.boundary.fun='tent';
options.boundary.choose=1:5;
options.boundary.color='red';

options.title=['Boundarys and Eigenfunctions of Tent Map with ',param.basis,' Basis (n=',num2str(param.n),')'];
% options.title=['Eigenfunctions of Tent Map with ',param.basis,' Basis (n=',num2str(param.n),',noise=',num2str(D),')'];
options.subp=[3,3];

options.save.enabled=true;
options.save.path='./temp';
options.save.pre=['Tent_eigen_',param.basis];
options.save.suffix='.png';
end

%% Logistic map
function [fun,param,options]=Logistic_map(options)
alpha=4;
D=0;
fun=@(x)awgn(alpha.*x.*(1-x),10*log10(1/D));

param.dim=1;
param.n=1000;
param.m=[1,2,3,4,5,6,7,8,9];%4;%[2,3,4,8,10,16,20,50,100];
param.phase=[0,1];
param.x0=linspace(param.phase(1),param.phase(2),param.n);
param.times=1;
param.basis='natural';

param.natural.enabled=true;%自然基
param.natural.x0=rand;

options.multim.enabled=true;%支持不同的基函数数量
options.multim.deal='real';
options.multim.choose=1;

options.boundary.enabled=true;
options.boundary.fun='logistic';
options.boundary.choose=1:5;
options.boundary.color='red';

options.title=['Boundarys and Eigenfunctions of Logistic Map with ',param.basis,' Basis (n=',num2str(param.n),')'];
%options.title=['Eigenfunctions of Logistic Map with ',param.basis,' Basis (n=',num2str(param.n),',noise=',num2str(D),')'];
options.subp=[3,3];

options.save.enabled=true;
options.save.path='./temp';
options.save.pre=['Logistic_eigen_',param.basis];
options.save.suffix='.png';
end

%% Intermittency Map
function [fun,param,options]=Intermittency(options)
fun=@(x)x./(1-x).*(x<=0.5)+(1-x)./x.*(x>0.5);

param.dim=1;
param.n=1000;
param.m=4;%[1,2,3,4,5,6,7,8,9];%[2,3,4,8,10,16,20,50,100];
param.phase=[0,1];
param.x0=linspace(param.phase(1),param.phase(2),param.n);
param.times=1;
param.basis='natural';%'natural'

param.natural.enabled=true;%自然基
param.natural.x0=rand;

options.multim.enabled=false;%支持不同的基函数数量
options.multim.deal='real';
options.multim.choose=1;

options.title=['Eigenfunctions of Intermittency Map with ',param.basis,' Basis (n=',num2str(param.n),')'];
% options.title=['Eigenfunctions of Tent Map with ',param.basis,' Basis (n=',num2str(param.n),',noise=',num2str(D),')'];
options.subp=[2,2];

options.save.enabled=false;
options.save.path='./temp';
options.save.pre=['Tent_eigen_',param.basis];
options.save.suffix='.png';
end

%% Default options
function options=default_options
options.subp=[3,3];

options.save.enabled=false;
options.save.path='./temp';
options.save.suffix='.png';
end