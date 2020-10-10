close all
fontsize=20;
issave=true;


%subplot(221)
figure
rect_draw(fontsize)
if issave
    saveas(gcf,'./temp/function_basis_rect.png')
end

%subplot(222)
figure
gauss_draw(fontsize)
if issave
    saveas(gcf,'./temp/function_basis_gauss.png')
end

%subplot(223)
figure
fourier_draw(fontsize)
if issave
    saveas(gcf,'./temp/function_basis_four.png')
end

%subplot(224)
figure
legendre_draw(fontsize)
if issave
    saveas(gcf,'./temp/function_basis_legen.png')
end
%set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);

function rect_draw(fontsize)
n=1000;m=4;
x0=linspace(0,1,n);
textpos=[linspace(1/2/m,1-1/2/m,m);1.8*ones(1,m)];
for i=1:m
    hold on
    g=Rectangle(i,m);
    plot(x0,g(x0));
    text(textpos(1,i)-0.03,textpos(2,i),['g_{',num2str(i),'}(x)'],'FontSize',fontsize)
end
%title('Rectangle Function Basis (m=4)')
sciformat;
end

function gauss_draw(fontsize)
n=1000;m=6;
x0=linspace(0,1,n);
textpos=[linspace(1/2/m,1-1/2/m,m);1.8*ones(1,m)];
for i=1:m
    hold on
    g=Gauss(i,m);
    plot(x0,g(x0));
    text(textpos(1,i)-0.03,textpos(2,i),['g_{',num2str(i),'}(x)'],'FontSize',fontsize)
end
%title('Gauss Function Basis (m=6)')
set(gca,'FontSize',fontsize);
sciformat;
end

function fourier_draw(fontsize)
n=1000;m=5;
x0=linspace(0,1,n);
textpos=[0.3,0.7;
    0.447,-1.338;
    0.325,1.259;
    0.519,1.372;
    0.688,0.985];
texts{1}='\sqrt{2}/2';
texts{2}='\sqrt{2}\cdot\cos{2\pi x}';
texts{3}='\sqrt{2}\cdot\sin{2\pi x}';
texts{4}='\sqrt{2}\cdot\cos{2\pi 2x}';
texts{5}='\sqrt{2}\cdot\sin{2\pi 2x}';
%textpos=[linspace(1/2/m,1-1/2/m,m);1.8*ones(1,m)];
color=lines(4);
color1=color(1,:);
color2=color(2,:);
color3=color(4,:);
colors=[color1;color2;color2;color3;color3];
for i=1:m
    hold on
    g=Fourier(i);
    plot(x0,g(x0),'Color',colors(i,:));
    %text(textpos(i,1)-0.04,textpos(i,2)+0.01,['$g_{',num2str(i),'}(x)=',texts{i},'$'],'interpreter','latex');
    text(textpos(i,1)-0.04,textpos(i,2)+0.01,['$g_{',num2str(i),'}(x)$'],'interpreter','latex','FontSize',fontsize);
end
%title('Fourier Function Basis (m=5)')
sciformat;
end

function legendre_draw(fontsize)
n=1000;m=5;
x0=linspace(-1,1,n);
% textpos=[
%     -1,0.707;
%     -0.8659,-1.06;
%     -0.904,1.147;
%     -0.8859,-0.7657;
%     -0.994,1.996;
%     -0.7578,0.9824];
textpos=[
    -1,0.707;
    -0.8659,-1.06;
    -0.904,1.147;
    -0.4635,0.835;
    -0.007,0.7951;
    0.3033,0.8085];
for i=0:m
    hold on
    g=Legendre(i);
    plot(x0,g(x0));
    text(textpos(i+1,1),textpos(i+1,2),['P_{',num2str(i),'}(x)'],'FontSize',fontsize)
end
%title('Legendre Function Basis (m=6)')
sciformat;
end

function g = Rectangle(i,m)
% x is a number
if(i==m)
    g=@(x)(x>=(i-1)/m & x<=i/m)*sqrt(m); %x坐标，第i个基函数
else
    g=@(x)(x>=(i-1)/m & x<i/m)*sqrt(m);
end
end



function g= Gauss(i,m)
dj=1/m/2;
xj=linspace(1/2/m,1-1/2/m,m);
c=1/sqrt(dj*sqrt(pi));
g=@(x,m)c*exp(-(x-xj(i)).^2./(2.*dj.^2));
end

function g= Fourier(i)
if i==1
    g=@(x)ones(1,length(x))*sqrt(2)/2;
elseif mod(i,2)==0
    g=@(x)sqrt(2)*cos(2*pi*floor(i/2).*x);
else
    g=@(x)sqrt(2)*sin(2*pi*floor(i/2).*x);
end
end

function g= Legendre(i)
getFirst=@(x)x(1,:);
g=@(x)getFirst(legendre(i,x,'norm'));
end