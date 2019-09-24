function B=Henon_Boundary(A,iter,iterinv)
if nargin==1
    iter=1;
    iterinv=1;
elseif nargin==2
    iterinv=iter;
end

a=1.4;
b=0.3;
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
f_inv=@(x,y)y/b;
g_inv=@(x,y)x-1+a/b/b*y.*y;
% syms x y;[X,Y]=solve(y+1-1.4*x*x==x,0.3*x==y);
% x1=0.6313544770895047116815602338357;
% x2=-1.1313544770895047116815602338357;
% y1=0.18940634312685141350446807015071;
% y2=-0.33940634312685141350446807015071;
% 
% Attr=load('./data/Henon_attractors_data_xy.mat'); % 吸引子数据载入
% x_attr=Attr.x;y_attr=Attr.y;
% Peri=load('./data/Henon_period_orbrits_P_1_0.3_-1_-0.3.mat');
% P=Peri.P;

B=A;
for i=1:length(A(:,1))
    x=A(i,1);y=A(i,2);
    x_temp=x;y_temp=y;
    for j=1:iter
        x_temp=f(x_temp,y_temp);y_temp=g(x_temp,y_temp);
        B=[B;x_temp,y_temp];
    end
    x_temp=x;y_temp=y;
    for j=1:iterinv
        x_temp=f_inv(x_temp,y_temp);y_temp=g_inv(x_temp,y_temp);
        B=[B;x_temp,y_temp];
    end
end