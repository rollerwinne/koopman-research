% clear;
% load('data0.mat');
% for i=2:9
%     for j=1:i-1
%         for k=1:length(X{i})
%             for l=1:length(X{j})
%                 if abs(X{i}(k)-X{j}(l))<1e-10
%                     X{i}(k)=NaN;
%                 end
%             end
%         end
%     end
% end
% for i=1:9
%     index=find(isnan(X{i}));
%     X{i}(index)=[];
% end
% % X{3}(2)=[];

x_k=x(1,:);
x_l=x(end,:);
for m=1:20
x0=linspace(1/2/m,1-1/2/m,m);

K=Logistic_G_Gauss(x_k,x0,m);
clf
figure
set(gcf,'outerposition',get(0,'screensize'));
for j=1:m
    hold on
    hhh=plot(x_k,K(:,j));
end
title(['Gauss Basis Function of m=',num2str(m)])
str=['.\temp\Logistic_Gauss_basis_m',num2str(m)];
saveas(hhh,[str,'.fig'])
saveas(hhh,[str,'.png'])
end
function G= Logistic_G_Gauss(x,xj,m)
dj=1/m/2;
G=zeros(length(x),m);
for i=1:length(x)
    for j=1:m
        G(i,j)=exp(-(x(i)-xj(j))^2/(2*dj^2));
    end
end
end