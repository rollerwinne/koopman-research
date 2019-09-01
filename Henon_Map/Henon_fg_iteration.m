% clear;clc
% x0=x0(i),y0=y0(i),T=T,Dx=k*v(1),Dy=k*v(2),times=dot,delta=0.5;
% save('.\data\temp.mat','x0','y0','T','Dx','Dy','times','delta')
%  load('.\data\temp.mat')
function [X,Y]=Henon_fg_iteration(x0,y0,T,Dx,Dy,times,delta)
dist=@(A,B)((A(1)-B(1))^2+(A(2)-B(2))^2)^0.5;
t=1e-5;
X=[x0];
Y=[y0];
for i=2:times
    count=0;
    while true
        x1=x0+t*Dx;
        y1=y0+t*Dy;
        [x2,y2]=Henon_fg_n(x1,y1,T*(i-1),1);
%         for j=1:T*(i-1)
%             x_temp=f_inv(x1,y1);
%             y_temp=g_inv(x1,y1);
%             x1=x_temp;
%             y1=y_temp;
%         end
        d=dist([X(i-1),Y(i-1)],[x2,y2]);
        if d<delta
            break;
        else
            t=t*0.1;
            count=count+1;
        end
        disp(['i=' num2str(i) '; count=' num2str(count) '; t=' num2str(t)]);
    end
    X=[X,x2];
    Y=[Y,y2];
end
end


