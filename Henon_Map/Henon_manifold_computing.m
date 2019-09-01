%function [XX,YY]=Henon_manifold_computing(x0,y0,Dx,Dy,T,t,n,dot)
%
% x0=x0(i),y0=y0(i),Dx=k*v(1),Dy=k*v(2),T=T,t=t,n=5,dot=2;
% save('.\data\temp.mat','x0','y0','Dx','Dy','T','t','n','dot');
% load('.\data\temp.mat');
%
distance=@(A,B)((A(1)-B(1))^2+(A(2)-B(2))^2)^0.5;
insert=@(A,n,a)[A(1:n),a,A(n+1:end)];
delete=@(A,n)[A(1:n-1),A(n+1:end)];
x1=x0+t*Dx;
y1=y0+t*Dy;
[x2,y2]=Henon_fg_n(x1,y1,1,1);
X{1}(1:n)=linspace(x1,x2,n);
Y{1}(1:n)=linspace(y1,y2,n);
delta_x=(x2-x1)/(n-1);
delta_y=(y2-y1)/(n-1);
delta=(delta_x^2+delta_y^2)^0.5;

for i=1:dot
    for j=1:length(X{i})
        [x_temp,y_temp]=Henon_fg_n(X{i}(j),Y{i}(j),T,1);
        X_temp(j)=x_temp;
        Y_temp(j)=y_temp;
    end
    %dist_01=ones(1,length(X_temp)-1); %判断距离二值化矩阵
    while true %sum(dist_01)>0 %使得所有距离都满足条件
        
        for j=1:length(X_temp)-1 %计算距离矩阵与二值化距离矩阵
            dist(j)=distance([X_temp(j),Y_temp(j)],[X_temp(j+1),Y_temp(j+1)]); %相邻点距离
            if dist(j)>delta
                dist_01(j)=1; %大于设定距离的置为1
            elseif dist(j)<=delta
                dist_01(j)=0; %小于设定距离的置为0
            end
        end
        
        if sum(dist_01)==0 %如果都满足设定距离 则循环结束
            break;
        end
        
        for j=find(dist_01==1) %1:length(X_temp)-1 %大于设定距离的,插值处理
            %if dist_01(j)==1
                X_aver=mean( [X{i}(j),X{i}(j+1)] );
                Y_aver=mean( [Y{i}(j),Y{i}(j+1)] );
                [X_aver_f,Y_aver_f]=Henon_fg_n(X_aver,Y_aver,T,1);
                X{i}=insert(X{i},j,X_aver);
                Y{i}=insert(Y{i},j,Y_aver);
                X_temp=insert(X_temp,j,X_aver_f);
                Y_temp=insert(Y_temp,j,Y_aver_f);
            %end
        end
        
    end %此循环结束后 X_temp 与 Y_temp 都满足距离要求
    X{i+1}=X_temp;
    Y{i+1}=Y_temp;
end

% 最后把所有的X{i}串到一起
XX=X{1};
YY=Y{1};
for i=2:dot+1
    XX=[XX,X{i}(2:end)];
    YY=[YY,Y{i}(2:end)];
end
%end
