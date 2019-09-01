clear;clc;
n=1000;part=4;
x0=linspace(0,1,n);
f1=@(x)(x<=0.5).*2.*x+(x>0.5).*(2-2*x);
f2=@(x)4*x.*(1-x);

for j=1:part
    temp1=round((j-1)/part*n)+1;
    temp2=round(j/part*n);
    x0_part=x0(temp1:temp2);
    x1=f1(x0_part);
    x2=f2(x0_part);
    C1=zeros(1,part);C2=C1;
    for i=1:length(x1)
        which=whichpart(part,x1(i));
        C1(which)=C1(which)+1;
        which=whichpart(part,x2(i));
        C2(which)=C2(which)+1;
    end
    D1(j,:)=C1/(n/part);
    D2(j,:)=C2/(n/part);
end

partstr=(1:part)/part;
X_node=cellstr(num2str(partstr'));
figure(1)
%set(gcf,'outerposition',get(0,'screensize'));
G1=digraph(D1,X_node);
hh=plot(G1,'EdgeLabel',G1.Edges.Weight);
axis square
saveas(hh,['.\temp\Logistic_graph_theory_part',num2str(part),'_degree1.png'])
figure(2)
%set(gcf,'outerposition',get(0,'screensize'));
G2=digraph(D2,X_node);
hh=plot(G2,'EdgeLabel',G2.Edges.Weight);
axis square
saveas(hh,['.\temp\Logistic_graph_theory_part',num2str(part),'_degree2.png'])

disp('The matrix is:');
D1
[f1,d1]=eig(D1)
D2
[f2,d2]=eig(D2)


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