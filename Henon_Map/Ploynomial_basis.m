function G=Ploynomial_basis(x,power)
%x=rand(4,1);dimention=3;power=2;
dimention=length(x(1,:));
n=length(x(:,1));
G=[];
for i=0:power
    if i==0
        G=[G,ones(n,1)];
    else
        P=Polynomial_fun(dimention,i,1);
        G_temp=[];
        for j=1:n
            G_temp=[G_temp;P(x(j,:))];
        end
        G=[G,G_temp];
    end
end
end