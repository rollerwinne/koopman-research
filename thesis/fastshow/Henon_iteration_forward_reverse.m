function Henon_iteration_forward_reverse
a=1.4;
b=0.3;
f=@(x,y)y+1-a.*x.*x;
g=@(x,y)b*x;
f_inv=@(x,y)y/b;
g_inv=@(x,y)x-1+a/b/b*y.*y;
[x,y]=meshgrid(linspace(-1.5,1.5,100));
for m=20
    %subplot(3,3,(m+1)/3)
    hold on
    for i=1:m
        x_temp=f_inv(x,y);
        y_temp=g_inv(x,y);
        x=x_temp;
        y=y_temp;
    end
    for i=1:m
        x_temp=f(x,y);
        y_temp=g(x,y);
        x=x_temp;
        y=y_temp;
    end
    plot3(x,y,0.04*ones(100,100),'b.')
    axis([-1.5 1.5 -1.5 1.5])
    %title(m)
end
% suptitle('Iteration of Henon Map')
% fullscreen;
% saveas(gcf,'./temp/Henon_iter_reverse_forward.png')
end