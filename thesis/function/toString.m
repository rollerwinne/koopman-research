function str=toString(m,split)
str='';
for i=1:length(m)-1
    str=[str,num2str(m(i)),split];
end
str=[str,num2str(m(end))];
end