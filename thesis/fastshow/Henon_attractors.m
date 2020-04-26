clear
a=1.4;
b=0.3;
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
% [X,Y]=solve('y+1-1.4*x*x=x','0.3*x=y');
x1=0.6313544770895047116815602338357;
x2=-1.1313544770895047116815602338357;
y1=0.18940634312685141350446807015071;
y2=-0.33940634312685141350446807015071;
% [x,y]=meshgrid(-1.5:0.01:1.5);
[x,y]=meshgrid(linspace(-1.5,1.5,100));
x0=-1.4;y0=-1;
for i=1:30
    figure(i)
    hold on
    plot(x,y,'b.')
    plot(x0,y0,'r*')
    x_temp=f(x0,y0);
    y_temp=g(x0,y0);
    x0=x_temp;
    y0=y_temp;
    %axis([-1.5 1.5 -1.5 1.5])
    x_temp=f(x,y);
    y_temp=g(x,y);
    x=x_temp;
    y=y_temp;
end
x=x(:);
y=y(:);
idx=find(isinf(x) | isinf(y));
x(idx)=[];y(idx)=[];

save('./data/Henon_attractors_n6245.mat','x','y','x1','x2','y1','y2');