function [Z]=project4_period_F(X,Y,f,g)

%% 定义n维向量Z(x,y),
global t;
Z=zeros(2*t,1);
for i=1:t
    if(i==1)
    Z(2*i-1,:)=X(i,:)-f(X(t,:),Y(t,:));
    Z(2*i,:)=Y(i,:)-g(X(t,:),Y(t,:));
    else Z(2*i-1,:) = X(i,:)-f(X((i-1),:),Y((i-1),:));
         Z(2*i,:) = Y(i,:)-g(X((i-1),:),Y((i-1),:));
    end
end