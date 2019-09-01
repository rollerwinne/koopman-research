clear %;clf,close all
a=1.4;
b=0.3;
%% Henon map equation
f=@(x,y)y+1-a*x.*x;
g=@(x,y)b*x;
dfgxy=@(x,y)[-2*a*x,1;b,0];
%% Initial point
% x1=[-1,0.2];
% x2=[1,-0.2];
x1=[1,0.3];
x2=[-1,-0.3];
% x1=[1,1];
% x2=[-1,-1];
% x1=[-1,1];
% x2=[1,-1];
load(['.\data\Henon_period_orbrits_P_' num2str(x1(1)) '_' num2str(x1(2)) '_' num2str(x2(1)) '_' num2str(x2(2)) '.mat']);

T=10;fig=6;
subplot(2,3,fig)
count=size(P{1,T},1); 
if count==0;ylabel('null'); end
for choose=1:count
ylabel({['T=' num2str(T) ',choose=' num2str(choose) '/' num2str(count) '/' num2str(P{1}(T))]})
X=P{1,T}(choose,mod(1:end,2)==1);
Y=P{1,T}(choose,mod(1:end,2)==0);
Z=0.02*ones(1,length(X));
hold on
h=plot3(X,Y,Z,'o','color','black','LineWidth',1,'MarkerSize',10)
drawnow
pause(2)
set(h,'visible','off');
end