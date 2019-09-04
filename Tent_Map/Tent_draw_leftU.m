clear;clc;close all
n=1000;p=1;

setup.function='RectangleL';
setup.complex='real';
setup.bin=0;
setup.figurenum=10;

x=Tent_x(n,p);
%x=Tent_x_noise(n,p,0.001);
x0=linspace(0,1,n);
peak_num=[];
for m=100
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
    for i=1:min(setup.figurenum,length(h))
        [pks,locs,deep]=Tent_findpeaks_leftU(n,m,F,D,K,x_marker,h,setup,0,i,1);
%         [pks,locs,deep]=Logistic_findpeaks3_2(n,m,F,D,U,h,setup,0.00,i,2,0);
        peak_num=[peak_num;m,length(pks),deep];
    end
%     for i=1:min(setup.figurenum,length(h))
%         figure
%         set(gcf,'outerposition',get(0,'screensize'));
%         for j=1%:4
%             %subplot(2,2,j)
%             switch j
%                 case 1
%                     E=real(K*F2(:,h(i)));
%                     hh=plot(x0,E);
%                     ylabel('real');
%                 case 2
%                     E=imag(F2(:,h(i)));
%                     hh=stem(x0,E,'.');
%                     ylabel('imaginary');
%                 case 3
%                     E=abs(F2(:,h(i)));
%                     hh=stem(x0,E,'.');
%                     ylabel('abs');
%                 case 4
%                     E=angle(F2(:,h(i)));
%                     hh=stem(x0,E,'.');
%                     ylabel('angle');
%             end
%         end
%         d_abs=abs(D(h(i)));
%         d_angle=angle(D(h(i)))/pi*180;
%         str1=['n=',num2str(n),'; m=',num2str(m)];
%         str2=[num2str(d_abs) ' б╧' num2str(d_angle) 'бу'];
%         suptitle({str1;str2});
%         str=['.\temp\Logistic_',setup.function,'_',setup.complex,'_n',num2str(n),'m',num2str(m),'_figure',num2str(i)];
% %         saveas(hh,[str,'.png'])
% %         saveas(hh,[str,'.fig'])
%     end
end
% close all