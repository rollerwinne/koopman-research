%% clean workspace
clc;
clear;
close all;
%% plot data
fprintf('plot data... \n');
X = load('ex2x.dat');
y = load('ex2y.dat');
m = length(y);
figure;
plot(X,y,'o');
%% gradient decent
fprintf('Runing gradient decent... \n');
X = [ones(m,1),X];
theta_SGD = zeros(2, 1);
theta_GD = zeros(2, 1);
theta_SVRG = zeros(2, 1);

Iteration = 2000;
alpha = 0.015;
alpha1 = 0.025;

[theta ,J]= StochasticGD(X, y, theta_SGD, alpha, Iteration);
[theta1 ,J1]= GradinentDecent(X, y, theta_GD, alpha, Iteration);
[theta2 ,J2]= SVRG(X, y, theta_SVRG, alpha1);

fprintf('SGD: %f %f\n',theta(1),theta(2));
fprintf('GD: %f %f\n',theta1(1),theta1(2));
fprintf('SVRG: %f %f\n',theta2(1),theta2(2));

hold on;
plot(X(:, 2), X*theta, 'r-');
plot(X(:, 2), X*theta1, 'g-');
plot(X(:, 2), X*theta2, 'b-');
legend('','SGD','GD','SVRG');

x_j = 1:1:20;
figure;
hold on;
plot(x_j, J, 'b-');
plot(x_j, J1, 'g-');
plot(x_j, J2, 'r-');
legend('SGD','GD','SVRG');
xlabel('epoch')
ylabel('loss')