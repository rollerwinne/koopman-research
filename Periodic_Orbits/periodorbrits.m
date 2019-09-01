clear;clf,close all
% alpha=4;
% f=@(x)alpha.*x.*(1-x);
% df=@(x)alpha-2*alpha*x;
f=@(x)1-2*abs(x-1/2);
df=@(x)(x>=0 & x<0.5)*2+(x==0.5)*0+(x>0.5 & x<=1)*(-2);
for n=2:10
    X=[];
    [A,num]=symbolperiod(n);
    for i=1:num
        x=zeros(n,1);
        x=x+(0.25+0.5*A(i,:))';
        sum=0;
        while true
            F=[];
            DF=[];
            for j=2:n
                F=[F;x(j)-f(x(j-1))];
                DF=[DF;zeros(1,j-2),-df(x(j-1)),1,zeros(1,n-j)];
            end
            F=[x(1)-f(x(n));F];
            DF=[1,zeros(1,n-2),-df(x(n));DF];
            delta=-DF\F;
            sum=sum+1;
            if norm(delta,2)<1e-4 || sum>1000;
                break;
            else
                x=x+delta;
            end
        end
        if sum<=1000
            X=[X;x'];
        else
            disp('²»ÊÕÁ²');
        end
    end
    P{n}=X;
end

C=[];
for i=2:10
   subplot(3,3,i-1)
   count=size(P{1,i},1);
   C=[C;count];
   title(['T=' num2str(i) ',count=' num2str(count)])
   hold on
   s=jet(count);
   for j=1:count
        plot(P{1,i}(j,:),f(P{1,i}(j,:)),'*','color',s(j,:))
   end
end

% choose=10;
% count=size(P{1,choose},1);
% for j=1:count
%     H=[];
%     for i=1:9
%         subplot(3,3,i)
%         hold on
%         s=jet(count);
%         h=plot(P{1,choose}(j,:),0.00*(j-count/2),'*','color','red');%,s(j,:))
%         H=[H,h];
%     end
%     drawnow;
%     if j==count
%         break;
%     end
%     pause(0.5);
%     for i=1:9
%         delete(H(:,i));
%     end
% end