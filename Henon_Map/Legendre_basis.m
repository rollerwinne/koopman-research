function G=Legendre_basis_1d(x,power)
%x=rand(4,1);dimention=3;power=2;
[n,dimention]=size(x);
G=[];
for i=0:power
    if i==0
        G=[G,ones(n,1)];
    else
        P=legendre(i,x,'norm');
        G_temp=P(1,:)';
        G=[G,G_temp];
    end
end
end
% n=5;
% X=-1:0.01:1;
% s=jet(n);
% for i=1:n
%     P=legendre(i,X,'norm');
%     for j=1
%         hold on
%         plot(X,P(j,:),'color',s(i,:))
%     end
% end
% 
% trapz(X,P(1,:).*P(1,:))