function rz=inter(x,y,z,T,choose)
period=load('./data/Henon_period_orbrits_P_1_0.3_-1_-0.3.mat'); % 周期轨道数据载入
P=period.P;
xp=P{1,T}(choose,mod(1:end,2)==1);
yp=P{1,T}(choose,mod(1:end,2)==0);
rz=zeros(1,10);
for i=1:length(xp)
    zp=interp2(x,y,z,xp(i),yp(i));
    rz(i)=relative(z,zp);
end
end

function rz=relative(Z,zp)
rz=(zp-min(Z(:)))/(max(Z(:))-min(Z(:)))*100;
end