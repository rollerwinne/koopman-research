clear;clc
data=load('tent_boundary_x0.66.mat');
X=data.X;
l=length(X);
for i=l:-1:2
    num=length(X{i-1});
    X{i}(num+1:end)=[];
end

% for i=3:l
%     num=length(X{i});
%     X{i}(1:num/2)=[];
% end
save('tent_boundary_norepeat_x0.66.mat','X');