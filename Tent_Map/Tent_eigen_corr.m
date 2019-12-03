function c=Tent_eigen_corr(K1,K2,L1,L2,n,m1,m2)
x0=linspace(0,1,n);
[h1,h2,F1,F2,D1,D2]=eigen(K1,K2,L1,L2);
c=zeros(length(h1),length(h2));
for i=1:length(h1)
    for j=1:length(h2)
        A1=real(K1*F1(:,h1(i)));
        A2=real(K2*F2(:,h2(j)));
        c(i,j)=corr2(A1,A2);
    end
end
hc=c(:);
[~,idx]=sort(abs(hc));
if length(hc)<9
    st=1;
else
    st=length(hc)-8;
end
fc=zeros(length(hc),1);
idx(st:end);
fc(idx(st:end))=1;
fc=reshape(fc,length(h1),length(h2));
%fc=(abs(c)>=abs(th)-1e-15);
k=1;
figure;
for i=1:length(h1)
    for j=1:length(h2)
        if (k<10 && fc(i,j)~=0)
            subplot(3,3,k)
            k=k+1;
            A1=real(K1*F1(:,h1(i)));
            A2=real(K2*F2(:,h2(j)));
            plot(x0,A1)
            hold on
            plot(x0,A2)
            title(['corr=',num2str(c(i,j)),'; \lambda_1=',num2str(D1(h1(i))),'; \lambda_2=',num2str(D2(h2(j)))])
            %legend(['m=',num2str(m1)],['m=',num2str(m2)])
        end
    end
end
suptitle(['m=',num2str(m1),'; m=',num2str(m2)])
end

function [h1,h2,F1,F2,D1,D2]=eigen(K1,K2,L1,L2)
% K1=K_pre;
% K2=K;
% L1=L_pre;
% L2=L;
U1=pinv(K1)*L1;
U2=pinv(K2)*L2;
[F1,D1]=eig(U1);
[F2,D2]=eig(U2);
D1=diag(D1);
D2=diag(D2);
h1=find(abs(D1)>1e-5 & abs(D1)<1.3 & abs(imag(D1))<1e-13);
h2=find(abs(D2)>1e-5 & abs(D2)<1.3 & abs(imag(D2))<1e-13);
% h1=1:length(D1);
% h2=1:length(D2);
end