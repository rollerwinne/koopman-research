clear
f=@(x)sin(x(:,1)).*cos(x(:,2));
[X,Y]=meshgrid(linspace(-2*pi,2*pi,100));
XX=X(:);YY=Y(:);
Z=reshape(f([XX,YY]),100,100);
mesh(X,Y,Z)
xlabel('x'),ylabel('y')

clear
f=@(x)sin(x(1)).*cos(x(2));
x0=[1.3,1.4];
x=fmin_GD(f,x0,0.00001,10000) 