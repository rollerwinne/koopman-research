clear all;clc;close all
rho=28;sigma=10;beta=8/3;
%f=@(x,y,z)[sigma*(y-x);x.*(rho-z)-y;x.*y-beta*z];
f=@(t,x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
x0=[-1,3,4];
dt=0.01;
tspan=[0,100,dt];
[T,X]=rk_4(f,tspan,x0);
% plot3(X(:,1),X(:,2),X(:,3))
% view([1,-1,0])
% xlabel('x'),ylabel('y'),zlabel('z')
fixedpoint1=[sqrt(beta*(rho-1)),sqrt(beta*(rho-1)),rho-1];
fixedpoint2=[-sqrt(beta*(rho-1)),-sqrt(beta*(rho-1)),rho-1];

Distance=@(x,y)sqrt(sum((x-y).^2)); % 1 row *2
LineSpeed=@(x,dt)sqrt(sum((x(2,:)-x(1,:)).^2))/dt; % 2 row
AngularSpeed=@(x,x0,dt)LineSpeed(x,dt)/Distance(x(1,:),x0); % 2 row
Acceleration=@(x,dt)(LineSpeed(x(2:3,:),dt)-LineSpeed(x(1:2,:),dt))/dt;% 3 row
KineticEnergy=@(m,x,dt)m/2*LineSpeed(x,dt).^2; % 2 row
str={'LineSpeed';'AngularSpeed';'Acceleration';'KineticEnergy'};
for i=1:4
    subplot(2,2,i)
    for j=2:length(T)-1
        if i==1
            Q(i,j)=LineSpeed(X(j:j+1,:),dt);
        elseif i==2
            if Distance(X(j,:),fixedpoint1)>Distance(X(j,:),fixedpoint2);
                Q(i,j)=AngularSpeed(X(j:j+1,:),fixedpoint1,dt);
            else
                Q(i,j)=AngularSpeed(X(j:j+1,:),fixedpoint2,dt);
            end
        elseif i==3
            Q(i,j)=Acceleration(X(j-1:j+1,:),dt);
        elseif i==4
            Q(i,j)=KineticEnergy(1,X(j:j+1,:),dt);
        end
    end
    hh=scatter3(X(2:end-1,1),X(2:end-1,2),X(2:end-1,3),3,Q(i,2:length(T)-1)','filled');
    title(str{i})
    colorbar
    colormap(jet)
    view([1 -1 0])
    axis equal
end