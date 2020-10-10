clc;clear;close all
figure(1)
Henon_attractors_draw('blue',0);
s=hsv(4);
for i=1:4
    Henon_boundary_draw(i-1,0,0,s(i,:),false,true);
end
xlabel('x'),ylabel('y')
sciformat
set(gcf,'outerposition',[0,0,560,560])
savesci('Henon_boundary_forward')
%title('Forward Iteration of Henon Boundary')
%saveas(gcf,'./temp/Henon_boundary_forward.png');

figure(2)
Henon_attractors_draw('blue',0);
s=hsv(4);
for i=1:4
    Henon_boundary_draw(-i+1,0,0,s(i,:),false,true);
end
xlabel('x'),ylabel('y')
sciformat
set(gcf,'outerposition',[0,0,560,560])
savesci('Henon_boundary_reverse')
%title('Reverse Iteration of Henon Boundary')
%saveas(gcf,'./temp/Henon_boundary_reverse.png');

% figure(3)
% Henon_attractors_draw('blue',true);
% Henon_boundary_draw(0,0,0,'red');
% title('Boundary of Henon Map')
%saveas(gcf,'./temp/Henon_boundary.png');
