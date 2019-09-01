clear;close all;clc
alpha=4;n=1000000;part=100;
f=@(x)alpha.*x.*(1-x);
root=@(x)roots([-alpha,alpha,-x]);
x=linspace(0,1,n);

X=reshape(x,n/part,part);
X=X';
A=zeros(part);
for i=1:part
    for j=1:size(X,2)
        root_temp=root(X(i,j));
        s_temp=[whichpart(part,root_temp(1)),whichpart(part,root_temp(2))];
        t_temp=ones(1,2)*i;
        %         if s_temp(2)==1 && t_temp(2)==10
        %             error('dd')
        %         end
        A(s_temp(1),t_temp(1))=A(s_temp(1),t_temp(1))+1;
        A(s_temp(2),t_temp(2))=A(s_temp(2),t_temp(2))+1;
    end
end
A=A/(n/part)/2;

figure(1)
set(gcf,'outerposition',get(0,'screensize'));
partstr=(1:part)/part;
X_node=cellstr(num2str(partstr'));
G=digraph(A,X_node);
hh=plot(G,'EdgeLabel',G.Edges.Weight);
axis square
str11='.\Logistic_graph\';
str22=['n',num2str(n),'part',num2str(part)];
str33='\graph';
saveas(hh,[str11,str22,str33,'.png'])
% figure(2)
% G2=graph(upper(A),X_node);
% %G.Edges.NormWeight = G.Edges.Weight/sum(G.Edges.Weight);
% plot(G2,'EdgeLabel',G2.Edges.Weight)
% % plot(G)
% axis square

% figure(3)
% rho=1;
% theta=linspace(0,2*pi,part);
% poly_x=rho.*cos(theta);
% poly_y=rho.*sin(theta);
% gplot(A,[poly_x',poly_y'],'-*')
% axis square


[F,D]=eig(A);
D=diag(D);
[~,index]=sort(abs(D));
h=index(end-10:end);
for j=1:length(h)
    figure(j+1)
    set(gcf,'outerposition',get(0,'screensize'));
    for i=1:4
        switch i
            case 1
                H=real(F(:,h(j)));
            case 2
                H=imag(F(:,h(j)));
            case 3
                H=abs(F(:,h(j)));
            case 4
                H=angle(F(:,h(j)));
        end
        subplot(2,2,i)
        hh=stem((1:part)/part,H,'.');
        %hh=bar((1:part)/part,H);
    end
    d_abs=abs(D(h(j)));
    d_angle=angle(D(h(j)))/pi*180;
    str1=['n=',num2str(n),'; part=',num2str(part)];
    str2=[num2str(d_abs) ' ¡Ï' num2str(d_angle) '¡ã'];
    suptitle({str1;str2});
    str3=['\graph_eigenfunctions_figure',num2str(j)];
    saveas(hh,[str11,str22,str3,'.png'])
end


function which=whichpart(part,x)
for i=1:part
    if i==part
        if x>=(i-1)/part && x<=i/part
            which=i;
        end
    else
        if x>=(i-1)/part && x<i/part
            which=i;
        end
    end
end
end