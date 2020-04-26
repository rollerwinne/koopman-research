function fullscreen(h)
if nargin<1
    set(gcf,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
else
    set(h,'outerposition',get(0,'screensize')-[0,0,1440*0.3,900*0.2]);
end
end

