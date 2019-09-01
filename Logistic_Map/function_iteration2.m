clear;clc;close all
figure_num=5;
F1=@(x)0.5+0.5*cos(2*pi*x);
F2=@(x)(x<=0.25).*1+(x>0.75).*1;
F3=@(x)(x<=0.5).*1;
F4=@(x)(x>=0.5).*1;

f1=@(x)(x<=0.5).*2.*x+(x>0.5).*(2-2*x);
f2=@(x)4*x.*(1-x);
x=0:0.001:1;
F=F4;

figure
set(gcf,'outerposition',get(0,'screensize'));
for i=0:figure_num
    subplot(figure_num+1,2,2*(i)+1)
    hh=stem(x,F(ff(f1,x,i)),'.');
    subplot(figure_num+1,2,2*(i)+2)
    hh=stem(x,F(ff(f2,x,i)),'.');
end
str=['.\temp\function_iteration2_figure4'];
saveas(hh,[str,'.fig']);
saveas(hh,[str,'.png']);

function fff=ff(f,x,n)
x_temp=x;
for i=1:n
    x_temp=f(x_temp);
end
fff=x_temp;
end
% for i=1:figure_num
%     subplot(figure_num,2,2*(i-1)+1)
%     plot(x,y1);
%     subplot(figure_num,2,2*(i-1)+2)
%     plot(x,y2);
%     y1=f1(y1);
%     y2=f2(y2);
% end
