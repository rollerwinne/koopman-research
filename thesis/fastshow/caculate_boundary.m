clear;clc
p=9;syms x;
X{1}=[2/3];
for i=2:p
    X{i}=[];
    for j=1:length(X{i-1})
        y=tent_solve(X{i-1}(j));
        X{i}=[X{i},y];
    end
end

function x=tent_solve(m)
temp=(1-m)/2;
x1=1/2-temp;
x2=1/2+temp;
x=[x1,x2];
end

function x=logistic_solve(m)
temp=sqrt(1-m)/2;
x1=1/2-temp;
x2=1/2+temp;
x=[x1,x2];
end