function [ theta_old, J_history ] = SVRG( X, y, theta, alpha )
theta_old = theta;
n = length(y);
J_history = zeros(20,1);
m = 2 * n;
for i = 1:20
     theta_ = theta_old;
     Mu = 1/n *  X' * (X*theta_ - y);
     theta_0 = theta_;
     for j = 1:m
         index = randi(n);
         GD_one = (X(index, :) * theta_0 - y(index)) * X(index, :)';
         GD_ = (X(index, :) * theta_ - y(index)) * X(index, :)';
         theta_t = theta_0 - alpha * (GD_one - GD_ + Mu);
         theta_0 = theta_t;
     end
     J_history(i) = ComputeCost(X, y, theta_t);
     theta_old = theta_t;
end
end