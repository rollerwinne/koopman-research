idx=findpeaks2(X(:),Y(:),Z(:));
plot3(X(:),Y(:),Z(:),'b.')
hold on;plot3(X(idx),Y(idx),Z(idx),'ro')

function idx=findpeaks2(x,y,z,err)
if nargin<4
    err=1e-1;
end
idx=[];
for i=1:length(x)
    set=[];
    for j=1:length(x)
        if distance(x(i),x(j),y(i),y(j))<err
            set=[set,i];
        end
    end
    set
    for j=set
        if z(i)==min(z(set)) || z(i)==max(z(set))
            idx=[idx,i];
        end
    end
end
end

function d=distance(x1,y1,x2,y2)
d=sqrt((x1-x2)^2+(y1-y2)^2);
end