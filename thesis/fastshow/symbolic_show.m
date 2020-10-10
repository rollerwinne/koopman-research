clear;close all
n=1000;
f=@(x)4.*x.*(1-x);
x0=linspace(0,1,n);
plot(x0,f(x0),'k');
xlabel('x_n'),ylabel('x_{n+1}');
hold on

syms x
x1=double(vpasolve(4*x*(1-x)-0.5,x));
plot([0,1],[0,1],'b--');
plot([0.5,0.5],[0,1],'r--');
plot([x1(1),x1(1)],[0,f(x1(1))],'r--');
plot([x1(2),x1(2)],[0,f(x1(2))],'r--');
plot([x1(1),x1(2)],[f(x1(1)),f(x1(2))],'r--');

str={'00','01','11','10'};
x_pos=[0.08,0.32,0.68,0.92];
for i=1:4
    text(x_pos(i),0.05,['''',str{i},''''],'FontSize',20,'HorizontalAlignment','center')
end
sciformat
set(gcf,'outerposition',[440   378   560   560])
%saveas(gcf,'./temp/Logistic_symbolic.png');