close all

path = '../rk4_cooling/Out3.txt';
path2 = '../rk4_cooling/macro_laguerre/Out4.txt';
[matrix, nPoints, nCols] = ReadTxt(path,false);
[matrix2, nPoints2, nCols2] = ReadTxt(path2,false);
x = matrix(:,2);
y = matrix(:,3);
z = matrix(:,4);
x2 = matrix2(:,2);
y2 = matrix2(:,3);
z2 = matrix2(:,4);

set(gcf, 'Position', [50,100,1280,600]) %window pos
tiledlayout(1,3);
c1 = nexttile;
plot3(x,y,z,''), grid on, title Mariana, xlabel x, ylabel y, zlabel z
axis([-1 8 -15 15 -15 15])
c2 = nexttile;
plot3(x2,y2,z2,'c'), grid on, title Manuela, xlabel x, ylabel y, zlabel z
axis([-1 8 -15 15 -15 15])
c3 = nexttile;
[X,Y] = meshgrid(-1:.1:1); 
Z=sin(X.*Y);
surf(X,Y,Z), xlabel x, ylabel y, zlabel z;
title('$$Sin(x\cdot y)$$','interpreter','latex');







function [matrix, nPoints, nCols] = ReadTxt(path, head)
    file = fopen(path, 'r');
    lts=0;
    if (head) lts=1; end
    matrix = dlmread(path,'',lts,0);
    nPoints = size(matrix,1);
    nCols = size(matrix,2);
    fclose(file);
end
