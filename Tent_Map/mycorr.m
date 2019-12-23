function XY=mycorr(x,y)
XY=0;X=0;Y=0;
xa=mean(x);
ya=mean(y);
for i=1:length(x)
    XY=XY+(x(i)-xa)*(y(i)-ya);
end
XY=XY/(length(x));
end