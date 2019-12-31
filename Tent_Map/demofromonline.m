clear; clc;
[X, Y, Z] = peaks(30);
position = [13,26; 20,24; 16,27; 16, 21; 16,24; 9,19; 9,15; 7,17; 12,17; 9,17; 14,16; 14,11; 
            11,12; 16,12; 14,13; 17,11; 16,5; 13,8; 20,8; 17,8; 22,19; 22,12; 20,16; 25,16; 22,16];
p1 = position(:, 1);
p2 = position(:, 2);
x = zeros(1,length(p1));
y = zeros(1, length(p1));
z = zeros(1, length(p1));
for i = 1 : length(p1)
    x(i) = X(p1(i), p2(i));
    y(i) = Y(p1(i), p2(i));
    z(i) = Z(p1(i), p2(i));
end
subplot(1,2,1)
surf(X, Y, Z)
hold on
plot3(x, y, z, 'or', 'MarkerFaceColor', 'r');
hold off
ylim([-2, 2]);
[Xi,Yi,Zi] = griddata(x,y,z,linspace(min(x),max(x), 30)',linspace(min(y),max(y), 30),'v4');
subplot(1,2,2)
surf(Xi,Yi,Zi)
hold on
plot3(x, y, z, 'or', 'MarkerFaceColor', 'r');
hold off