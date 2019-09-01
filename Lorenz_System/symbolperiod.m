function [A,num]=symbolperiod(n)
% 给出一个n,输出不同的符号周期轨道及符号周期轨道的个数
A=[];
for i=1:2^n-2
    A=[A;bitget(i,n:-1:1)];
end

row=1;
while row<size(A,1)
    a=A(row,:);
    for j=1:n
        a=circshift(a',1)';
        delete=[];
        for k=row+1:size(A,1)
            if a==A(k,:)
                delete=[delete,k];
            end
        end
        for l=length(delete):-1:1
            A(delete,:)=[];
        end
    end
    row=row+1;
end

b=1:n;
c=b(mod(n,b)==0);
c([1,end])=[];

row=1;
while row<=size(A,1)
    for i=1:length(c)
        times=n/c(i);
        flag=0;
        for j=2:times
            if sum(A(row,1:c(i))~=A(row,(j-1)*c(i)+1:j*c(i)))>0
                flag=1;
            end
        end
        if flag==0
            A(row,:)=[];
            break;
        end
    end
    row=row+1;
end
num=size(A,1);