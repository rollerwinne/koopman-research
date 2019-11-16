function [f,seq,sx]=Tents_function_low(n,D)
seq=linspace(0,1,n);
block=seq(2);
k=1/seq(2);
str='@(x)';awgnstr=str;
if (nargin>1)
    awgnstr=[awgnstr,'awgn('];
end
for i=1:n-1
    if (i==n-1) %最后一个取双闭区间
        area=['.*(x>=',num2str(i-1),'*block&x<=',num2str(i),'*block)'];
    else
        area=['.*(x>=',num2str(i-1),'*block&x<',num2str(i),'*block)'];
    end
    if(mod(i,2)==1) %奇偶表达式不同
        if (i==n-1 || i==n-2)
            equ=['(-',num2str((i-1)/2),'+k*x/2)'];
        else
            equ=['(-',num2str(i-1),'+k*x)'];
        end
    elseif (mod(i,2)==0)
        if (i==n-1 || i==n-2)
            equ=['(',num2str(i/2),'-k*x/2)'];
        else
            equ=['(',num2str(i),'-k*x)'];
        end
        
    end
    str=[str,'+',equ,area];
    awgnstr=[awgnstr,'+',equ,area];
end
if(nargin>1)
    awgnstr=[awgnstr,',10*log10(1/D))'];
end
disp(awgnstr)
disp(str)
f=eval(awgnstr);
x0=seq(1:end-1)+block/2;
sx=[];
for i=2:n-1
    sx=[sx;fsolve(eval([str,'-',num2str(seq(i))]),x0)];
end
end