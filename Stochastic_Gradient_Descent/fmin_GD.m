function x=fmin_GD(f,x0,yita,iter)
x=x0;
    disp(x);
for i=1:iter
    for j=1:length(x0)
        dx(j)=DFx(f,x,j,0.00001);
    end
    x=x-yita*dx;

end
end

function dx=DFx(f,x,i,delta_x)
dx=(f(x+[zeros(i-1,1);delta_x;zeros(length(x)-i,1)])-f(x))/delta_x;
end