function [fun,l,f_split]=Polynomial_fun(d,p,flag) 
% d:dimention (d>=1)
% p;power (p>=1)
% flag=1(row) or flag=0(queue)
x_str=[];
x_addstr=[];
for i=1:d
    x_str=[x_str,' x',num2str(i)];
    if i==1
        x_addstr=[x_addstr,'x',num2str(i)];
    else
         x_addstr=[x_addstr,'+x',num2str(i)];
    end
end
str=['syms',x_str];
eval(str);
str=['f=expand((',x_addstr,')^',num2str(p),');'];
eval(str);
f_str=char(f);
f_str(isspace(f_str)) = [];
f_split=regexp(f_str,'+','split');
f_split=f_split';
f_funstr=[];
l=length(f_split);
for i=1:l
    if isstrprop(f_split{i}(1),'digit')
        f_where=strfind(f_split{i},'*');
        f_split{i}(1:f_where(1))=[];
    end
    f_split{i}=strrep(f_split{i},'*','.*');
    f_split{i}=strrep(f_split{i},'^','.^');
	[x_start,x_end]=regexp(f_split{i},'x[0-9]+');
    while ~isempty(x_start)
        f_split{i}=[f_split{i}(1:x_start(1)),'(',f_split{i}(x_start(1)+1:x_end(1)),')',f_split{i}(x_end(1)+1:end)];
        [x_start,x_end]=regexp(f_split{i},'x[0-9]+');
        %break;
    end
    if i==1
        f_funstr=f_split{i};
    elseif flag==1
        f_funstr=[f_funstr,',',f_split{i}]; % in a row
    else
        f_funstr=[f_funstr,';',f_split{i}]; % in a queue
    end
end
str=['fun=@(x)[',f_funstr,'];'];
eval(str);
end