function [ theta,J_history ] = StochasticGD( X, y, theta, alpha, num_iter )
m = length(y);
J_history = zeros(20, 1);
temp = 0;
n = 0;
for iter = 1:num_iter
     temp = temp + 1;
     index = randi(m);
     theta = theta -alpha *  (X(index, :) * theta - y(index)) * X(index, :)';
     if temp>=100
         temp = 0;
         n = n + 1;
         J_history(n) = ComputeCost(X, y, theta);
     end
end
end