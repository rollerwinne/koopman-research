function K=Henon_attractors_flow(n,iter)
a=1.4;b=0.3;
f=@(x,y)y+1-a*x.*x;              
g=@(x,y)b.*x;

x_fix=((b-1)+((b-1).^2+4*a).^(1/2))*((2*a).^-1);
y_fix=b*x_fix;
A=[-2*a*x_fix,1;b,0];
[V,~]=eig(A);     
k=V(2,1)/V(1,1);

x_step=linspace(-0.02,0.02,n);
x0=x_fix+x_step';
y0=k*(x0-x_fix)+y_fix;


for i=1:iter
x0_temp=f(x0,y0);
y0_temp=g(x0,y0);
x0=x0_temp;
y0=y0_temp;
end
K=zeros(2*n,1);
K(1:2:end)=x0;
K(2:2:end)=y0;

% figure;
% axis([-1.5,1.5,-0.4,0.4]);
% hold on;
% for i=1:length(x0)
%     plot(x0(i),y0(i),'b.');
%     drawnow;
% end
end