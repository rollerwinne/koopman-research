clearvars -except F D x_function_num;clc;close all
tic;timestart=char(datetime('now'));
dj=5;%m=4;
choose='complex';  % 'real' or 'complex'
bin=0;
tspan=[0,100,0.01];
x0=[-1,3,4];
[t,x]=Lorenz_Koopman_x(tspan,x0);
x_k=x(1:end-1,:);
x_l=x(2:end,:);
%[F,D,U,x_function_num]=Lorenz_Koopman_U(x_k,x_l,dj);
%save('.\data\Lorenz_Matrix_data_FD_t0_100_0.02dj5.mat','F','D','x_function_num');
%load('.\data\Lorenz_Matrix_data_FD_t0_100_0.01dj5.mat');
if bin==1
    F(F>=0)=1;F(F<0)=-1; % binarization
end
if strcmp(choose,'real')==1
    h=find(abs(D)>0.9 & abs(D)<1.1 & imag(D)==0); % find real eigenvalues
elseif strcmp(choose,'complex')==1
    h=find(abs(D)>0.9965 & abs(D)<1.0035 & imag(D)>-1e-6 ); % find complex eigenvalues
end

F_start=100;
for i=1:min(6,length(h))
    figure(2)
    set(gcf,'outerposition',get(0,'screensize'));
    subplot(6,1,i)
    
    d_abs=abs(D(h(i)));
    d_angle=angle(D(h(i)))/pi*180;
    str1=['n=',num2str(length(x_k(1:end,1))),'; m=',num2str(x_function_num(1)),'*',num2str(x_function_num(2)),'*',num2str(x_function_num(3)),...
        '=',num2str(x_function_num(1)*x_function_num(2)*x_function_num(3)),'; d_{j}=',num2str(dj),''];
    str2={[num2str(d_abs)];['б╧',num2str(d_angle),'бу']};
    
    F_X=F_start:length(x_k);
    F_temp=abs(F(F_X,h(i)));
    F_findnum=floor(length(F_temp)*0.1);
    [~,F_index]=sort(F_temp);
    F_which=F_index(end-F_findnum:end);
    
    hold on
    h2=stem(t(F_X),F_temp,'.');
    h2=plot(t(F_X(F_which)),F_temp(F_which),'r*');
    suptitle({str1})
    ylabel(str2)
    
    figure(1)
    set(gcf,'outerposition',get(0,'screensize'));
    subplot(2,3,i)
    hold on
    XX=x_k(F_X,1);YY=x_k(F_X,2);ZZ=x_k(F_X,3);
    %plot3(XX,YY,ZZ)
    h1=scatter3(XX,YY,ZZ,3,F_temp,'filled');
    colorbar
    colormap(jet)
    view([1 -1 0])
    axis equal
    h1=plot3(XX(F_which),YY(F_which),ZZ(F_which),'r*');%,'color','blue','Markersize',3);
    str1=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
    title(str1)
    
end
str3=['.\temp\Lorenz_eigenfunctions1d_findmax',choose,'_figure1.png'];
str4=['.\temp\Lorenz_eigenfunctions1d_findmax',choose,'_figure2.png'];
saveas(h1,str3);
saveas(h2,str4);
