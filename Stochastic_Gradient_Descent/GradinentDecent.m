function [ theta, J_history ] = GradinentDecent( X, y, theta, alpha, num_iter )
m = length(y);
J_history = zeros(20, 1);
i = 0;
temp = 0;
for iter = 1:num_iter
     temp = temp +1;
     theta = theta - alpha / m * X' * (X*theta - y);
     if temp>=100
         temp = 0;
         i = i + 1;
         J_history(i) = ComputeCost(X, y, theta);
     end
end
end