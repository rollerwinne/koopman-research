function G=Legendre_basis_3d(x,power)
%x=x_k;
%x_boundary=[-20,20;-30,30;0,50];
width=[20,30,25];
center=[0,0,25];
for i=1:3
    x_norm(:,i)=(x(:,i)-center(i))/width(i);
end
G=[];
for i=0:power
    if i==0
        P=legendre(0,x_norm(:,1),'norm')/sqrt(width(1))/sqrt(width(2))/sqrt(width(3));
        G=[G,P(1,:)'];
    else
        for nx=0:i
            for ny=0:i-nx
                PX=legendre(nx,x_norm(:,1),'norm')/sqrt(width(1));
                PY=legendre(ny,x_norm(:,2),'norm')/sqrt(width(2));
                PZ=legendre(i-nx-ny,x_norm(:,3),'norm')/sqrt(width(3));
                G_temp=PX(1,:)'.*PY(1,:)'.*PZ(1,:)';
                G=[G,G_temp];
            end
        end
    end
end
end