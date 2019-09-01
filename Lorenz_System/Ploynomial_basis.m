function [G,G_str,G_fun]=Ploynomial_basis(x,power)
%x=rand(6,3);power=3;
[n,dimention]=size(x);
G=[];
for i=0:power
    if i==0
        G=[G,ones(n,1)];
        G_str='1';
    else
        [P,~,split]=Polynomial_fun(dimention,i,1);
        G_str=[G_str;split];
        G_temp=[];
        for j=1:n
            G_temp=[G_temp;P(x(j,:))];
        end
        G=[G,G_temp];
    end
end
G_fun=funstr2handle(G_str);
end

function handle=funstr2handle(str)
for i=1:length(str)
    if i==1
        funstr=str{i};
    else
        funstr=[funstr,';',str{i}]; % in a queue
    end
end
str1=['@(x)[',funstr,'];'];
handle=eval(str1);
end