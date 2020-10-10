function savesci(str,h)
if nargin<2
    h=gcf;
end
saveas(h,['./temp/',str,'.png']);
%saveas(h,['./temp/',str,'.fig']);
%saveas(h,['./temp/',str,'.eps'],'psc2');
%print('-depsc',['./temp/',str,'.eps']);
end

