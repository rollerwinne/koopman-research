function [XX,YY]=Henon_manifold_computing4(x0,y0,Dx,Dy,T,t,n,dot)
%
% x0=x0(i),y0=y0(i),Dx=k*v(1),Dy=k*v(2),T=T,t=t,n=n,dot=dot;
% save('.\data\temp.mat','x0','y0','Dx','Dy','T','t','n','dot');
% clear
% load('.\data\temp.mat');
%
distance=@(A,B)((A(1)-B(1))^2+(A(2)-B(2))^2)^0.5;
insert=@(A,n,a)[A(1:n),a,A(n+1:end)];
u0_overbar_x=@(u0,u1,u2)real(u1+(u1-u2)/abs(u1-u2)*abs(u1-u0)); % use complex input
u0_overbar_y=@(u0,u1,u2)imag(u1+(u1-u2)/abs(u1-u2)*abs(u1-u0)); % use complex input
alpha=@(u0,u0_overbar,u1)2*asin(abs(u0_overbar-u0)/2/abs(u1-u0)); % use complex input
% delete=@(A,n)[A(1:n-1),A(n+1:end)];
x1=x0+t*Dx;
y1=y0+t*Dy;
[x2,y2]=Henon_fg_n(x1,y1,T,1);
X{1}(1:n)=linspace(x1,x2,n);
Y{1}(1:n)=linspace(y1,y2,n);
delta=distance([x1,y1],[x2,y2])/(n-1);
% delta_x=(x2-x1)/(n-1);
% delta_y=(y2-y1)/(n-1);
% delta=(delta_x^2+delta_y^2)^0.5;