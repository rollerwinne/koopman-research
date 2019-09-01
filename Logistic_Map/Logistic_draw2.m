clear;clc;close all
alpha=4;n=1000;p=1;

setup.function='Rectangle';
setup.complex='complex';
setup.bin=0;
setup.figurenum=9;

x=Logistic_x(alpha,n,p);
%x=Logistic_x_noise(alpha,n,p,0.01);
x0=linspace(0,1,n);
for m=10
    switch setup.function
        case {'GaussS'} % Sparse Gauss basic function
            %m=3;
            [F,D] = Logistic_U_Gauss_sparse(x,m);
        case {'Gauss'} % Gauss basic function
            m=1000;md=1000;
            [F,D] = Logistic_U_Gauss(x,m,md);
        case {'Fourier'} % Fourier basic function
            %m=100;
            [F,D] = Logistic_U_Fourier(x,m);
        case {'Rectangle'} % Rectangle basic function
            m=100;
            [F,D] = Logistic_U_Rectangle(x,m);
        case {'RectangleN'} % Rectangle natural basic function
            m=1000;
            [F,D] = Logistic_U_Rectangle(x,m,1); % if nargin>2 use natural basis
        case {'Legendre'}
            % To be finished
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
        h=find(abs(D)>0.7 & abs(D)<1.3 & imag(D)==0); % find real eigenvalues
    elseif strcmp(setup.complex,'complex')==1
        h=find(abs(D)>0.2 & abs(D)<1.8 & imag(D)>-1e-6 ); % find complex eigenvalues
    end
    
    for i=1:min(setup.figurenum,length(h))
        figure
        set(gcf,'outerposition',get(0,'screensize'));
        for j=1:4
            subplot(2,2,j)
            switch j
                case 1
                    E=real(F2(:,h(i)));
                    hh=stem(x0,E,'.');
                    ylabel('real');
                case 2
                    E=imag(F2(:,h(i)));
                    hh=stem(x0,E,'.');
                    ylabel('imaginary');
                case 3
                    E=abs(F2(:,h(i)));
                    hh=stem(x0,E,'.');
                    ylabel('abs');
                case 4
                    E=angle(F2(:,h(i)));
                    hh=stem(x0,E,'.');
                    ylabel('angle');
            end
        end
        d_abs=abs(D(h(i)));
        d_angle=angle(D(h(i)))/pi*180;
        str1=['n=',num2str(n),'; m=',num2str(m)];
        str2=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
        suptitle({str1;str2});
        str=['.\temp\Logistic_',setup.function,'_',setup.complex,'_n',num2str(n),'m',num2str(m),'_figure',num2str(i)];
%         saveas(hh,[str,'.png'])
%         saveas(hh,[str,'.fig'])
    end
end