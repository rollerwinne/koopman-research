function [J] = project4_period_FD(X)
%% 定义F(x,y)的n阶导
global t;
global k;  
b=1;
d=1;
J=eye(2*t);
for i=1:(t-1)
   x=X(i,:);
   aa=1+k*cos(2*pi*x);
   cc=k*cos(2*pi*x);
  J(2*i+1:2*i+2,2*i-1:2*i)=-[aa b;cc d];
end
   x=X(t,:);
   aa=1+k*cos(2*pi*x);
   cc=k*cos(2*pi*x);
  J(1:2,2*t-1:2*t)=-[aa b;cc d];
% ff=@(x,y)(x+y+k/(2*pi)*sin(2*pi*x));
% gg=@(x,y)(y+k/(2*pi)*sin(2*pi*x));
% syms xx
% syms yy
% a=diff(ff(xx,yy),xx);
% b=diff(ff(xx,yy),yy);
% c=diff(gg(xx,yy),xx);
% d=diff(gg(xx,yy),yy);
% J=eye(2*n);
% for i=1:(n-1)
%     aa=double(subs(a,[xx,yy],[X(i,1),Y(i,1)]));
%     bb=double(subs(b,[xx,yy],[X(i,1),Y(i,1)]));
%     cc=double(subs(c,[xx,yy],[X(i,1),Y(i,1)]));
%     dd=double(subs(d,[xx,yy],[X(i,1),Y(i,1)]));
%     J(2*i+1:2*i+2,2*i-1:2*i)=-[aa bb;cc dd];
% end
%     aa=double(subs(a,[xx,yy],[X(n,1),Y(n,1)]));
%     bb=double(subs(b,[xx,yy],[X(n,1),Y(n,1)]));
%     cc=double(subs(c,[xx,yy],[X(n,1),Y(n,1)]));
%     dd=double(subs(d,[xx,yy],[X(n,1),Y(n,1)]));
%     J(1:2,2*n-1:2*n)=-[aa bb;cc dd];
end

