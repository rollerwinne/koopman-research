function Gauss_function(x,m,md)

x0=linspace(0,1,m);
G=Logistic_G_Gauss(x,x0,m,md);
figure(1)
for i=1:m
    hold on
    plot(x,G(:,i))
end
title(['n=',num2str(length(x)),' m=',num2str(m),' dj=1/',num2str(md)])
G_sum=sum(G,2);
figure(2)
plot(x,G_sum)
ylim([0,max(G_sum)])
title(['n=',num2str(length(x)),' m=',num2str(m),' dj=1/',num2str(md)])
end

function G= Logistic_G_Gauss(x,xj,m,md)
dj=1/md;
G=zeros(length(x),m);
for i=1:length(x)
    for j=1:m
        G(i,j)=exp(-(x(i)-xj(j))^2/(2*dj^2));
    end
end
end

