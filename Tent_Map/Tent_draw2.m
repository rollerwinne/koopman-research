clear;clc;close all
n=1000;p=1;

setup.function='Rectangle';
setup.complex='real';
setup.bin=0;
setup.figurenum=10;

x=Tent_x(n,p);
x0=linspace(0,1,n);
peak_num=[];
for m=10
    switch setup.function
        case {'GaussS'} % Sparse Gauss basic function
            %m=3;
            [F,D,U] = Tent_U_Gauss(x,m);
        case {'Gauss'} % Gauss basic function
            m=1000;md=1000;
            [F,D] = Logistic_U_Gauss(x,m,md);
        case {'Fourier'} % Fourier basic function
            %m=100;
            [F,D] = Logistic_U_Fourier(x,m);
        case {'Rectangle'} % Rectangle basic function
            %m=100;
            [F,D,U] = Tent_U_Rectangle(x,m);
        case {'RectangleN'} % Rectangle natural basic function
            m=1000;
            [F,D] = Logistic_U_Rectangle(x,m,1); % if nargin>2 use natural basis
        case {'Legendre'}
            % To be finished
    end
    matstr=['.\data\Tent_data_FD_',setup.function,'_n',num2str(n),'_m',num2str(m),'.mat'];
    % save(matstr,'F','D');
    % load(matstr)
    
    if setup.bin==0
        F2=F;
    else
        F2(F>=0)=1;F2(F<0)=-1; % binarization
    end
    if strcmp(setup.complex,'real')==1
        h=find(abs(D)>0.00000001 & abs(D)<1.3 & imag(D)==0); % find real eigenvalues
    elseif strcmp(setup.complex,'complex')==1
        h=find(abs(D)>0.01 & abs(D)<1.3 & imag(D)>-1e-6 &imag(D)~=0); % find complex eigenvalues
    end
    
    for i=1:min(setup.figurenum,length(h))
        [pks,locs,deep]=Tent_findpeaks(n,m,F,D,h,setup,0,i,1);
        peak_num=[peak_num;m,length(pks),deep];
    end
end