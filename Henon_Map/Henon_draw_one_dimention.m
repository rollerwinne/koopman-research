clear;clc;%close all
tic;timestart=char(datetime('now'));
a=1.4;b=0.3;q=10000;n=1;nn=4;
m=50;md=45;
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
subplot(4,3,[4,7])
load('.\data\Henon_attractors_data_xy.mat'); % 吸引子数据载入
plot3(x,y,0.02*ones(1,length(x)),'.','color','yellow','Markersize',1)  %吸引子图像

x=0;
y=0;
X=x;
Y=y;
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
hh=plot3(X(1,1:8),Y(1,1:8),0.025*ones(1,length(X(1,1:8))),'.','color','black','Markersize',10);
%legend('Attractors','Trace:All','Trace:1-6','Location','Best','-')
xy_text={'1','2','3','4','5','6','7','8'};
text(X(1,1:8)+0.02,Y(1,1:8)+0.02,0.021*ones(1,length(X(1,1:8))),xy_text,'color','red')
view(0,90)
axis([-1.5 1.5 -1.5 1.5])
F_start=100;
for i=1:n
    [F,D,U] = Henon_U(X(i,F_start:q),X(i,F_start+1:q+1),Y(i,F_start:q),Y(i,F_start+1:q+1),m,md);
    save('.\data\Henon_Matrix_data1d_FD_n1e5m50md45.mat','F','D')
    h=find(abs(D)>0.8 & abs(D)<1.2 & imag(D)>-1e-6);
    for j=1:min(nn,length(h))
        subplot(nn,3,[3*j-1 3*j])
        d_abs=abs(D(h(j)));
        d_angle=angle(D(h(j)))/pi*180;
        hh=stem(F_start:q,imag(F(:,h(j))),'.');
        %hh=plot(1:q,F(:,h(j)));
        str1=[num2str(d_abs) ' ∠' num2str(d_angle) '°'];
        title(str1)
    end
end
str1=['.\temp\Henon_eigenfunctions1d_complex_figure.fig'];
str2=['.\temp\Henon_eigenfunctions1d_complex_figure.png'];
%if ~isempty(h)
saveas(hh,str1);
saveas(hh,str2);
%end

timeuse=toc;
timeend=char(datetime('now'));
filename=mfilename('fullpath');
[filepath,name,~]=fileparts(filename);
address=java.net.InetAddress.getLocalHost;
%IPaddress=char(address.getHostAddress);
IPaddress=char(address);
subject=['Program: ''',name,'.m'' has been finished'];
content={['Time used: ',num2str(timeuse),' seconds.'];
    ['Start time: ',timestart,'.'];
    ['Finish time: ',timeend,'.'];
    ['Address: ',IPaddress,'.'];
    ['This is an auto E-mail from MATLAB.']}
qqmail2me(subject,content)