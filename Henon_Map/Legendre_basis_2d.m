function G=Legendre_basis_2d(x,power)
%x=rand(4,1);dimention=3;power=2;
width=[1.5,1.5];
center=[0,0];
for i=1:2
    x_norm(:,i)=(x(:,i)-center(i))/width(i);
end
G=[];
for i=0:power
    if i==0
        P=legendre(0,x_norm(:,1),'norm')/sqrt(width(1))/sqrt(width(2));
        G=[G,P(1,:)'];
    else
        for nx=0:i
            PX=legendre(nx,x_norm(:,1),'norm')/sqrt(width(1));
            PY=legendre(i-nx,x_norm(:,2),'norm')/sqrt(width(2));
            G_temp=PX(1,:)'.*PY(1,:)';
            G=[G,G_temp];
        end
    end
end
end

