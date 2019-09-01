clc;clear;%close all
rho=28;sigma=10;beta=8/3;
f=@(x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
period='0001';T=length(period);
% x1=[13,17,27];
% x2=[-13,-17,27];
x0=[-1,3,4];z0=27;
a=[0,0,1];
x1=Lorenz_Poincare_next_point(x0,z0,1);
J0=eye(3);
X0=[x1(:)',J0(:)'];
X1=Lorenz_Jacobian_next_point(X0,z0,T-1);
x=x1';xyz=[];
for i=1:T-1
    x=[x;X1(i,1:3)'];
end
xyz=[X0;X1];
count=0;
while true
    F1=[];DF1=[];DF2=[];DF3=[];
    for i=2:T
        xyz_temp=Lorenz_Jacobian_next_point(xyz(i-1,:),z0,1);
        F1=[F1;[xyz(i,1:3)-xyz_temp(1:3)]'];
        DF1=[DF1;[zeros(3,3*(i-2)),-reshape(xyz_temp(end-8:end),3,3),eye(3),zeros(3,3*(T-i))]];
        DF2=[DF2;[zeros(3,1*(i-1)),-f(xyz_temp(1:3)),zeros(3,1*(T-i))]];
        DF3=[DF3;[zeros(1,3*(i-1)),a,zeros(1,3*(T-i))]];
    end
    xyz_temp=Lorenz_Jacobian_next_point(xyz(T,:),z0,1);
    F1=[[xyz(1,1:3)-xyz_temp(1:3)]';F1];
    F2=zeros(T,1);
    DF1=[eye(3),zeros(3,3*(T-2)),-reshape(xyz_temp(end-8:end),3,3);DF1];
    DF2=[-f(xyz_temp(1:3)),zeros(3,1*(T-1));DF2];
    DF3=[a,zeros(1,3*(T-1));DF3];
    DF4=zeros(T);
    
    F=[F1;F2];
    DF=[DF1,DF2;DF3,DF4];
    delta=-DF\F;
    count=count+1;
    
    if norm(delta,2)<1e-8 || count>100
        break;
    else
        disp(['T=',num2str(T),',period=[',num2str(period),'],count=',num2str(count)]);
        x=x+delta(1:3*T);
        xyz(:,1)=x(mod(1:end,3)==1);
        xyz(:,2)=x(mod(1:end,3)==2);
        xyz(:,3)=x(mod(1:end,3)==0);
    end
end
if count<=100
    X=[X;x'];
    disp(['T=' num2str(T),',period=',period,',ÊÕÁ²']);
else
    disp(['T=' num2str(T),',period=',period,',²»ÊÕÁ²']);
end

