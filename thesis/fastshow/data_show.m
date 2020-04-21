clc;
f=@(x)1-2*abs(x-1/2);
g=@(x)4.*x.*(1-x);
d=dir('./data');
N=3:length(d);
for i=N
    str=d(i).name;
    data=load(str);
    flag=0;
    try
        data.X;
    catch
        flag=1;
    end
    if flag==0
        disp([str,':::::X']);
        cell2str(data.X);
        clear data.X
    else
        disp([str,':::::P']);
        cell2str(data.P);
        clear data.P
    end
end

function cell2str(X)
for i=1:length(X)
    disp(toString(X{i},','));
end
end

function str=toString(m,split)
if isempty(m)
    str='null';
    return;
end
str='[';
for i=1:length(m)-1
    str=[str,num2str(m(i),'%.4f'),split];
end
str=[str,num2str(m(end)),']'];
end