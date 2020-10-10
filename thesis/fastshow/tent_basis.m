close all;clc;

myfigure
basis2
savesci('Tent_U2')
myfigure
basis4
savesci('Tent_U4')
myfigure
U2
savesci('Koopman_eig2')
myfigure
U4
savesci('Koopman_eig4')

function U2
A=[.5,.5;.5,.5];
[F,D]=eig(A);D=diag(D);
x0=linspace(0-eps,1+eps,1000);
for m=1:2
    subplot(1,2,3-m)
    plot([0,1],[0,0],'r');
    hold on
    g=@(x)F(1,m).*(x<.5)+F(2,m).*(x>=.5);
    plot(x0,g(x0)*sqrt(2),'b','LineWidth',3);
    d_abs=abs(D(m));
    d_angle=angle(D(m))/pi*180;
    str=[num2str(d_abs),'∠' num2str(d_angle),'°'];
    title(str);
    sciformat(20)
    axis square
    axis([0-eps,1+eps,-1,1])
end
end

function U4
A=[.5,.5,0,0;0,0,.5,.5;0,0,.5,.5;.5,.5,0,0];
[F,D]=eig(A');D=diag(D);
F=[1,1,1;
    1,1,-1;
    1,-1,1;
    1,-1,-1];
D=[1,0,0]
x0=linspace(0,1,1000);
for m=1:3
    subplot(2,2,m)
    plot([0,1],[0,0],'r');
    hold on
    g=@(x)F(1,m).*(x<.25&x>=0)+F(2,m).*(x>=.25&x<.5)+F(3,m).*(x<.75&x>=.5)+F(4,m).*(x>=.75);
    plot(x0,real(g(x0)),'b','LineWidth',3);
    d_abs=abs(D(m));
    d_angle=angle(D(m))/pi*180;
    str=[num2str(d_abs),' ∠' num2str(d_angle),'°'];
    title(str);
    sciformat(20)
    axis square
    axis([0,1,-1,1])
end
end

function basis4
tightsub(3,3,4,0.7)
fun_basis(4,[1,2,3,4],false)
title({'rectangular basis (m=4)';''})
x1=[0.28,0.37];y1=[0.5,0.88];
for i=1:4
    textt([0.04+(i-1)*0.25,.9],['$g_',num2str(i),'(x)$'],15);
    annotation(gcf,'arrow',x1,y1+[0,-0.25]*(i-1),'LineWidth',1);
end


str={'''1000''','''0100''','''0010''','''0001'''};
for i=1:4
    tightsub(4,3,2+(i-1)*3,0.6)
    fun_basis(4,i,true)
    textt([-.5,.9],['$g_',num2str(i),'(x)$'],15)
    title(str{i})
    textt([1.5,0.5],'$\stackrel{U}{\longrightarrow}$',40)
end

num=[1,8;2,7;3,6;4,5];
for i=1:4
    tightsub(4,3,3+(i-1)*3,0.6)
    fun_basis(8,num(i,:),true)
end
end

function basis2
tightsub(3,3,4,0.7)
fun_basis(2,[1,2],false)
title({'rectangular basis (m=2)';''})
x1=[0.27,0.34];y1=[0.5,0.75];
for i=1:2
    textt([0.18+(i-1)*0.5,.9],['$g_',num2str(i),'(x)$'],15)
    annotation(gcf,'arrow',x1,y1+[0,-0.5]*(i-1),'LineWidth',1);
end

str={'''10''','''01'''};
for i=1:2
    tightsub(2,3,2+(i-1)*3,0.6)
    fun_basis(2,i,true)
    textt([-.4,.9],['$g_',num2str(i),'(x)$'],20)
    title(str{i})
    textt([1.17,0.5],'$\stackrel{U}{\longrightarrow}$',30)
end

num=[1,4;2,3];
for i=1:2
    tightsub(2,3,3+(i-1)*3,0.6)
    fun_basis(4,num(i,:),true)
end
end

function fun_basis(M,m,flag)
hold on;
x0=linspace(0-eps,1+eps,1000);
if flag
    i=1;
    for temp_m=m
        g{i}=Rectangle(temp_m,M);
        i=i+1;
    end
    sum=zeros(1,length(x0));
    for j=1:i-1
        sum=sum+g{j}(x0);
    end
    plot(x0,sum,'b','LineWidth',3);
else
    for temp_m=m
        g=Rectangle(temp_m,M);
        plot(x0,g(x0),'LineWidth',3);
    end
end

sciformat(20)
axis equal
axis([0,1,0,1])
end

function g = Rectangle(i,m)
if(i==m)
    g=@(x)(x>=(i-1)/m & x<=i/m)*sqrt(1);
else
    g=@(x)(x>=(i-1)/m & x<i/m)*sqrt(1);
end
end

function textt(x,str,fontsize)
text(x(1),x(2),str,'FontSize',fontsize,'Interpreter','Latex')
end

function textmatrix(x)
start_str='\begin{pmatrix}';
str='';
for i=1:length(x)-1
    for j=1:length(x)
        str=[str,num2str(x(i,j)),' & '];
    end
    str=[str,num2str(x(i,end)),' \\ '];
end
end_str='\end{pmatrix}';
matrix=[start_str,str,end_str];
textt([0,0],matrix,30);
end

function f = latex2(A, precision)
% 没输入精度参数时，默认精度为小数后4位
if nargin == 1
    precision = '4';
else
    precision = int2str(precision);
end
% 定义单一元素输出格式
out_num = [' %0.' precision 'f &'];
% 用作整数输出判断
z = zeros(1, str2num(precision) + 1);
z(1) = '.';
z(2 : end) = '0';
z = char(z);
% 求矩阵大小
[r c] = size(A);
nc = zeros(1, c);
nc(:) = 99;  % 存放character c
% 生成第一句Latex语句
out = sprintf('\\left(\n\t\\begin{array}{%s}', char(nc));
% 二重循环，用作生成整个矩阵的Latex语句
for i = 1 : r
    out = [out sprintf('\n\t')]; % 换行
    for j = 1 : c
        temp = sprintf(out_num, A(i, j));
        % 小数位皆为零时，把数取整。如1.0001取为1
        dot_position = find(temp == '.');
        if temp(dot_position : end - 2) == z
            temp = temp(1 : dot_position - 1);
            temp = [temp ' &'];
            % 要取整时，如有负号，则必须丢掉
            if temp(2) == '-'
                temp = [temp(1) temp(3 : end)];
            end
        end
        out = [out temp];
    end
    % 丢掉最后的'&'号
    out = out(1 : end - 1);
    % 行末加上'\\'号
    out = [out '\\'];
end
% 加上最后一句结束代码
out = [out,sprintf('\n\t\\end{array}\n\\right)')];
f = out;
end