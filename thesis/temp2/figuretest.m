clear
x=1:10;
y=5*x.^2+x;
subplot(1,2,1);plot(x,y,'r');
set(gca,'FontSize',20);
subplot(1,2,2);plot(x,y);
set(gca,'FontSize',20);