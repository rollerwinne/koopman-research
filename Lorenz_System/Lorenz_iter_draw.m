clear;close all;
n=5000;m=20;
rho=28;sigma=10;beta=8/3;
%f=@(x,y,z)[sigma*(y-x);x.*(rho-z)-y;x.*y-beta*z];
f=@(t,x)[sigma*(x(2)-x(1));
    x(1).*(rho-x(3))-x(2);
    x(1).*x(2)-beta*x(3)];
x0=[-1,3,4];
tf=0.01;t=100;
[T,X]=rk_4(f,[0,t,tf],x0);
%plot3(X(:,1),X(:,2),X(:,3));
x1=X(end-m-n+1:end,:);
x2=reshape(x1',(m+n)*3,1);
K=[];L=[];
for i=1:m
    K(:,i)=x2(3*i-2:3*(n+i-1));
    L(:,i)=x2(3*(i+1)-2:3*(n+i));
end
U=pinv(K)*L;
[F,D]=eig(U);
D=diag(D);
%% choose real or complex
choose='complex';
if strcmp(choose,'real')==1
    h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.8 & imag(D)>-1e-6 ); % find complex eigenvalues
end
if length(D)<=9
    h=1:length(D);
end
%%
figure_num=10;
for i=1:min(figure_num,length(h))
    figure(i)
    phi=K*F(:,h(i));
    for j=1:6
        switch(j)
            case 1
                Z=real(phi(1:3:end));
                str='x-real';
            case 2
                Z=abs(phi(1:3:end));
                str='x-abs';
            case 3
                Z=real(phi(2:3:end));
                str='y-real';
            case 4
                Z=abs(phi(2:3:end));
                str='y-abs';
            case 5
                Z=real(phi(3:3:end));
                str='z-real';
            case 6
                Z=abs(phi(3:3:end));
                str='z-abs';
            otherwise
                disp('error');
        end
        subplot(3,2,j)
        hh=scatter3(x1(1:n,1),x1(1:n,2),x1(1:n,3),3,Z,'filled');
        zlabel(str);
        colorbar
        colormap(jet)
        view([1 -1 0])
        axis equal
    end
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    str1=['n=',num2str(n),'; m=',num2str(m)];
    str2=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
    suptitle({str1;str2});
end