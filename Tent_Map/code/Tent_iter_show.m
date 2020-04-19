clear; close all;
digits(1024);
f=@(x)1-(2-1e-6)*abs(x-1/2);
n=60;
init=sqrt(3)/2;
hello=zeros(1,n);
for i=1:n
    hello(i)=init;
    init=f(init);
end
%hello=K(:,1);
Hello=f(hello);
s=jet(length(hello));
hold on
axis([0,1,0,1])
for i=1:length(hello)
    plot(hello(i),Hello(i),'*','Color',s(i,:))
    text(hello(i)+0.01,Hello(i)+0.01,num2str(i));
    if i>1
        plot(hello(i-1:i),Hello(i-1:i),'Color',[0.6,0.6,0.6])
    end
    pause(0.05)
    drawnow;
end
colorbar
colormap(jet)