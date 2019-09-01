clearvars -except F D;clc;%close all
tic;timestart=char(datetime('now'));
a=1.4;b=0.3;q=10000;n=1;nn=4;
m=50;md=45;
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
figure
subplot(4,3,[4,7])
load('.\data\Henon_attractors_data_xy.mat'); % 吸引子数据载入
plot3(x,y,0.02*ones(1,length(x)),'.','color','yellow','Markersize',1)  %吸引子图像

x=0;X=x;
y=0;Y=y;
title(['x1=',num2str(x),'; y1=',num2str(y)])
for i=1:q
    x_temp=f(x,y);
    y_temp=g(x,y);
    x=x_temp;
    y=y_temp;
    X=[X,x];
    Y=[Y,y];
end
hold on
plot3(X(1,:),Y(1,:),0.021*ones(1,length(X(1,:))),'.','color','blue','Markersize',3);
view(0,90)
axis([-1.5 1.5 -1.5 1.5])

F_start=100;
for i=1:n
    %[F,D,U] = Henon_U(X(i,F_start:q),X(i,F_start+1:q+1),Y(i,F_start:q),Y(i,F_start+1:q+1),m,md);
    %save('.\data\Henon_Matrix_data1d_FD_n1e5m50md45.mat','F','D')
    %load('.\data\Henon_Matrix_data1d_FD_n1e5m50md45.mat')
    h=find(real(D)>0 & abs(D)>0.85 & abs(D)<1.2 & imag(D)>-1e-6);
    for j=1:min(nn,length(h))
        subplot(nn,3,[3*j-1 3*j])
        hold on
        d_abs=abs(D(h(j)));
        d_angle=angle(D(h(j)))/pi*180;
        F_X=F_start:q;
        F_temp=abs(F(:,h(j)));
        F_findnum=floor(length(F_temp)*0.01);
        [~,F_index]=sort(F_temp);
        F_which=F_index(end-F_findnum:end);
        hh=stem(F_X,F_temp,'.');
        plot(F_X(F_which),F_temp(F_which),'r*');
        %hh=plot(1:q,F(:,h(j)));
        str1=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
        title(str1)
        if j==4
            subplot(4,3,[4,7])
            hold on
            XX=X(1,F_X);YY=Y(1,F_X);
            plot3(XX(F_which),YY(F_which),0.022*ones(1,length(XX(F_which))),'r*');%,'color','blue','Markersize',3);
            view(0,90)
            axis([-1.5 1.5 -1.5 1.5])
        end
    end
end
str1=['.\temp\Henon_eigenfunctions1d_findmax_complex_figure.fig'];
str2=['.\temp\Henon_eigenfunctions1d_findmax_complex_figure.png'];
%if ~isempty(h)
saveas(hh,str1);
saveas(hh,str2);
%end

%% send an E-mail to me
%[a,b]=qqmail2me(timestart,mfilename('fullpath'),[]); %程序开始时间、文件名、附件