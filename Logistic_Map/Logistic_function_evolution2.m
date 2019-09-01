clear;close all;clc
alpha=4;n=20000;part=200;p=1000;percent=5;
rate=1;%[0.5;0.5];
lamda=1;
f=@(x)alpha.*x.*(1-x);
root=@(x)roots([-alpha,alpha,-x]);
x=linspace(0,1,n);
g=@(x)gaussmf(x,[0.01,0]);
X=x;G=g(x);
G=[ones(1,n*percent/100),zeros(1,n*(100-percent)/100)];
figure(1)
stem(X,G,'.');
figure(2)
for j=1:p
    X_temp=[];G_temp=[];
    for i=1:length(X)
        if abs(G(i))<1e-16
            continue;
        else
            xroot=root(X(i));
            X_temp=[X_temp;xroot(randperm(2,1))];
            G_temp=[G_temp;rate*G(i)*lamda];
        end
    end
    X=X_temp;
    G=G_temp;
%     if j>4 && j<11
%         subplot(2,3,j-4)
%         stem(X,G,'.')
%         title(j)
%     end
end
stem(X,G,'.')
hold on
plot(0:1/part:1,G(1),'r*')
hold off
GG=[];

for i=1:part
    if i==part%length(x)
        which=find(X>=(i-1)/part & X<=i/part);
        w(i)=length(which);
        GG(i)=sum(G(which));
        GG2(i)=length(G(which));
    else
        which=find(X>=(i-1)/part & X<i/part);
        w(i)=length(which);
        GG(i)=sum(G(which));
        GG2(i)=length(G(which));
    end
end
figure(3)
set(gcf,'outerposition',get(0,'screensize'));
for i=1
    switch i
        case 1
            H=real(GG);
            H2=real(GG2);
        case 2
            H=imag(GG);
            H2=imag(GG2);
        case 3
            H=abs(GG);
            H2=abs(GG2);
        case 4
            H=angle(GG);
            H2=angle(GG2);
    end
    %subplot(2,2,i)
    %hh=fill([1:part,fliplr(1:part)],[zeros(size(H)),fliplr(H)],'k');
    %hh=yyaxis(linspace(0,1,part),H,linspace(0,1,part),H2);
    %yyaxis left
    %stem(linspace(0,1,part),H,'.');
    %yyaxis right
    hh=stem(linspace(0,1,part),H2,'.');%,'Visible','off');
    hold on
    plot(0:1/part:1,0,'r*')
    %xlim([0,1])
end
title(['n=',num2str(n),'; part=',num2str(part),'; p=',num2str(p)])
saveas(hh,['.\temp\','n',num2str(n),'part',num2str(part),'percent',num2str(percent),'%p',num2str(p),'_5.png'])
