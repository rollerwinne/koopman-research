%clear
%set=[];
%err=0.2;
%set=belong(2.1,set,err);
function auto_level(f,xp,yp,A,err)
bound=0.5;
res=[];
for i=1:length(xp)
    temp=xp(i);
    count=1;
    while(~equals(temp,bound,err))
        if count>15
            hh=text(xp(i),yp(i),'');
            break;
        end
        count=count+1;
        temp=f(temp);
    end
    if count<=15
        res(i)=count;
        hh=text(xp(i),yp(i),num2str(count));
    end
end
end


function auto_level2(f,x,err)
setx={};
sety={};
for i=1:length(x)
    temp=f(x(i));
    for j=1:length(x)
        if ( equals(x(j),temp,err))
            setx=add(x(i),setx,err);
            sety=add(x(j),sety,err);
        end
    end
end
arrx=cell2arr(setx);
arry=cell2arr(sety);
disp(arrx);
hold on
plot(arrx,0,'ro');
hold on
plot(arry,0.01,'ko');
end

% function auto_level2(f,x,err)
% for i=1:length(x)
%     temp1=f(x(i));
% end
% end

function res=equals(x,y,err)
if(abs(x-y)<err)
    res=1;
else
    res=0;
end
end

function set=add(x,set,err)
for i=1:length(set)
    for j=1:length(set{i})
        if (equals(x,set{i}(j),err))
            set{i}(end+1)=x;
            return;
        end
    end
end
set{end+1}=[x];
end

function res=belong(x,arr,err)
for i=1:length(arr)
    if (equals(x,arr(i))==1)
        res=i;
        return;
    end
end
res=-1;
end

function arr=cell2arr(cell)
arr=zeros(1,length(cell));
for i=1:length(cell)
    arr(i)=mean(cell{i});
end
end