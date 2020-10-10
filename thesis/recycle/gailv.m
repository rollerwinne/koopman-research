clear;close all
f=@(n,k)1./(1-(k-1)./n);%放回抽样
g=@(n,k)1+(k-1)./n;%不放回抽样 %有问题，这样算不对
n=10000;
for k=1:10000
    x(k)=k;
    y1(k)=f(n,k);
    y2(k)=g(n,k);
    y1_sum(k)=sum(y1);
    y2_sum(k)=sum(y2);
    %disp([num2str(k),':',num2str(f(10000,k))])
end
myfigure
subplot(221)
plot(x,y1,'r.','LineWidth',2)
title('放回抽样')
xlabel('抽签次数'),ylabel('单次抽奖次数的数学期望')
sciformat(15)

subplot(222)
plot(x,y2,'r.','LineWidth',2)
title('不放回抽样')
xlabel('抽签次数'),ylabel('单次抽奖次数的数学期望')
sciformat(15)

subplot(223)
plot(x,y1_sum,'r.','LineWidth',2)
title('放回抽样')
xlabel('抽签次数'),ylabel('总抽奖次数的数学期望')
sciformat(15)

subplot(224)
plot(x,y2_sum,'r.','LineWidth',2)
title('不放回抽样')
xlabel('抽签次数'),ylabel('总抽奖次数的数学期望')
sciformat(15)

suptitle(['抽签总数:',num2str(n),' 中奖人数:',num2str(k)])

% myfigure
% plot(x,y1,'r.','LineWidth',2)
% title('放回抽样')
% xlabel('抽签次数'),ylabel('单次抽奖次数的数学期望')
% sciformat
% savesci('y1')
% 
% myfigure
% plot(x,y2,'r.','LineWidth',2)
% title('不放回抽样')
% xlabel('抽签次数'),ylabel('单次抽奖次数的数学期望')
% sciformat
% savesci('y2')
% 
% myfigure
% plot(x,y1_sum,'r.','LineWidth',2)
% title('放回抽样')
% xlabel('抽签次数'),ylabel('总抽奖次数的数学期望')
% sciformat
% savesci('y1sum')
% 
% myfigure
% %subplot(224)
% plot(x,y2_sum,'r.','LineWidth',2)
% title('不放回抽样')
% xlabel('抽签次数'),ylabel('总抽奖次数的数学期望')
% sciformat
% savesci('y2sum')

%suptitle(['抽签总数:',num2str(n),' 中奖人数:',num2str(k)])