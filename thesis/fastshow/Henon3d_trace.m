clear;
A=-1.86;
C=0.03;
fun1=@(x)[x(2);
    x(3);
    0.72*x(1)+C*x(2)+A*x(3)-1.45*x(3).^2+0.515.*x(2).*x(3)-x(2).^2];
x0=[0.1,0,0];
%x0=rand(1,3);
n=10000;
x=zeros(n,3);
for i=1:1000
    x(i,:)=x0;
    x0=fun1(x0);
end
figure;
plot3(x(:,1),x(:,2),x(:,3),'.')

clear;
A=0.82;
C=2.06;
fun2=@(x)[x(2);
    x(3);
    0.5.*x(1)+A*x(3)+C*x(2)-2.25*x(3).^3-2*x(2).^3];
x0=[0.1,0,0];
%x0=rand(1,3);
n=10000;
x=zeros(n,3);
for i=1:1000
    x(i,:)=x0;
    x0=fun2(x0);
end
figure;
plot3(x(:,1),x(:,2),x(:,3),'.')