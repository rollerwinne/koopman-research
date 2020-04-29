function D1_boundary_draw(fun,Choose,yt,color)
X=getBoundary(fun);
hold on;
for choose=Choose
    plot(X{choose},yt,'*','Color',color)
end
end

function x=getBoundary(fun)
switch fun
    case {'tent'}
        X=load('tent_boundary_norepeat_x0.mat');
    case {'logistic'}
        X=load('logistic_boundary_norepeat_x0.mat');
    otherwise
        error('fun not support');
end
x=X.X;
end
