clear;clc;close all
n=1000;p=1;

setup.function='GaussL';
setup.complex='real';
setup.bin=0;
setup.figurenum=10;

x=Tent_x(n,p);
%x=Tent_x_noise(n,p,0.001);
x0=linspace(0,1,n);
peak_num=[];
for m=2
    switch setup.function
        case {'GaussL'}
            [F,D,U,K,L,x_marker] = Tent_U_Gauss_leftU(x,m);
        case {'RectangleL'}
            [F,D,U,K,L,x_marker] = Tent_U_Rectangle_leftU(x,m);
            rank(F)
    end
    matstr=['.\data\Logistic_data_FD_',setup.function,'_n',num2str(n),'_m',num2str(m),'.mat'];
    % save(matstr,'F','D');
    % load(matstr)
    
    if setup.bin==0
        F2=F;
    else
        F2(F>=0)=1;F2(F<0)=-1; % binarization
    end
    if strcmp(setup.complex,'real')==1
%         h=find(abs(D)>0.001 & abs(D)<1.3 & imag(D)==0); % find real eigenvalues
        h=1:length(D);
    elseif strcmp(setup.complex,'complex')==1
        %h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6 &imag(D)~=0); % find complex eigenvalues
        h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6); % find complex eigenvalues
    end
    figure;
    for i=1:min(setup.figurenum,length(h))
        [pks,locs,deep]=Tent_findpeaks_leftU(n,m,F,D,K,x_marker,h,setup,i);
%         [pks,locs,deep]=Logistic_findpeaks3_2(n,m,F,D,U,h,setup,0.00,i,2,0);
        peak_num=[peak_num;m,length(pks),deep];
    end
end
% close all