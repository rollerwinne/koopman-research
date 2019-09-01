function [tout,yout]=rk_4(odefile,tspan,y0)
ts=tspan;
t0=ts(1);
tf=ts(2);
yout=[];
tout=[];
y0=y0(:);
yout=[yout;y0'];
tout=[tout;t0];
if length(ts)==3
    h=ts(3);
else
    h=(ts(2)-ts(1))/100;
    tf=ts(2);
end
for t=[t0:h:tf-h]
    k1=h*feval(odefile,t,y0);
    k2=h*feval(odefile,t+h/2,y0+0.5*k1);
    k3=h*feval(odefile,t+h/2,y0+0.5*k2);
    k4=h*feval(odefile,t+h,y0+k3);
    y0=y0+(k1+2*k2+2*k3+k4)/6;
    yout=[yout;y0'];
    tout=[tout;t+h];
end

