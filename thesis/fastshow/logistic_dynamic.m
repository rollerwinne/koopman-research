clear;clc;close all
hold on
a=2;
A1=[0,a];B1=[2,a];C1=[1,a];
A2=[0,0];B2=[4,0];C2=[2,0];
A3=[0,-a];B3=[0,-a];C3=[2,-a];

plott(A1,B1);plott(A2,B2);plott(A3,C3);

quiverr(A1,A2);quiverr(B1,B2);quiverr(C1,C2);
quiverr(A2,A3);quiverr(B2,B3);quiverr(C2,C3)

textt(A1,'0');textt(B1,'1');textt(C1,'0.5');
textt(A2,'0');textt(B2,'2');textt(C2,'1');
textt(A3,'0');textt(B3,'0');textt(C3,'1');

texttt(A1,'A');texttt(B1,'B');texttt(C1,'C');
texttt(A2,'A''');texttt(B2,'B''');texttt(C2,'C''');
texttt(A3,'A''''');texttt(B3-[0.11,0]*2,'B''''');texttt(C3,'C''''');

eqstr1='T(x)=(T2\cdot T1)(x)=4x(1-x)';eqx1=[2.3,-a/3*2];
eqstr2='T_1(x)=2x';eqx2=[1-0.6,0+1];
eqstr3='T_2(x)=1-(x-1)^2';eqx3=[1-0.6,-a+1];
textb(eqx1,eqstr1);textb(eqx2,eqstr2);textb(eqx3,eqstr3);

axis off
set(gcf,'outerposition',[440,378,560,560])
sciformat
axis equal
savesci('Logistic_dynamic');

function plott(x1,x2)
plot([x1(1),x2(1)],[x1(2),x2(2)],'k-o','MarkerSize',5,'MarkerFaceColor','k','LineWidth',2.5);
end

function quiverr(x1,x2)
quiver(x1(1),x1(2),x2(1)-x1(1),x2(2)-x1(2),'Color','black','AutoScale',1,'AutoScale','off','MaxHeadSize',0.1,'LineWidth',.8);
end

function textb(x,str)
text(x(1),x(2),['$',str,'$'],'Color','black','FontSize',15,'Interpreter','latex');
end

function textt(x,str)
plot(x(1),x(2),'o','MarkerFaceColor','black','MarkerEdgeColor','black')
text(x(1)-0.12,x(2)-0.1,str,'Color','blue','FontSize',15)
end

function texttt(x,str)
text(x(1)+0.11,x(2)+0.11,str,'Color','red','FontSize',15,'HorizontalAlignment','center')
end

