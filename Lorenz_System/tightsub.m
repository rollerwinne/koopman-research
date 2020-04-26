function tightsub(M,N,i,weight)
w = 1.0/N; h = 1.0/M;
wi = weight*w; hi = weight*h;
[x,y] = ind2sub([N,N],i);
subplot('Position', [(x-1)*w+(1-weight)*w/2 1-y*h+(1-weight)*h/2 wi hi]);
end