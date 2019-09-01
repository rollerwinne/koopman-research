tic
A=[97.9748,59.4896,11.7418,8.5516,73.0331;
43.8870,26.2212,29.6676,26.2482,48.8609;
11.1119,60.2843,31.8778,80.1015,57.8525;
25.8065,71.1216,42.4167,2.9220,23.7284;
40.8720,22.1747,50.7858,92.8854,45.8849];
jordan(int16(A))
toc
A=U;
e=eig(A);
e(16:end)=0;
Tzz=unique(e(:));
Jor=[];
for i=1:length(Tzz)
    num=length(find(e==Tzz(i)));
    zhi=rank(eye(size(A))*Tzz(i)-A);
    [a,b]=size(A);
    x=a-zhi;
    js=num-x+1;
    for j=1:js
        [m,n]=size(Jor);
        Jor(m+1,n+1)=Tzz(i);
        if j>1
            Jor(m,n+1)=1;
        end
    end
    if x>1
        [m,n]=size(e);
        for k=1:(x-1)
            Jor(m+k,n+k)=Tzz(i);
        end
    end
end
Jor