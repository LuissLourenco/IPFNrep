close all 

a0=30;
r0=1;
%phi0=0;
p0=-10;
p=1;
l=3;
pasta = strcat('Out_a0_',num2str(a0),'_r0_',num2str(r0),'_p0_',num2str(p0),'_pl_',num2str(p),num2str(l),'/');
%pasta = strcat('Out_a0_',num2str(a0),'_phi0_',num2str(phi0),'_p0_',num2str(p0),'_pl_',num2str(p),num2str(l),'/');
%pasta = 'Outputs/';
Files=dir(pasta);
Names = {Files.name};
Names=Names(3:end);
Names = string(Names);
n = size(Names,2);
res = cell(n,1);

for k=1:n
    [matrix] = ReadTxt(strcat(pasta,Names(k)),false);
    res{k} = mat2cell(matrix, size(matrix,1), size(matrix,2));
end




%TRAJECTORIES
f=figure;
tiledlayout(1,3)
set(gcf,'Position',[50 50 1200 500])
ax1 = nexttile;
grid on, xlabel x, ylabel y, zlabel z;
view(ax1,3);
hold on;
for k=5
    matrix = cell2mat(res{k});
    x = matrix(:,2);
    y = matrix(:,3);    
    z = matrix(:,4);
    plot3(x,y,z)
end

ax2 = nexttile;
title(['a_0=',num2str(a0),'    r_0=',num2str(r0),'    p_0_x=',num2str(p0),'    pl=',num2str(p),num2str(l)]);
grid on, xlabel x, ylabel y, zlabel z;
view(ax2,3);
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    x = matrix(:,2);
    y = matrix(:,3);    
    z = matrix(:,4);
    plot3(x,y,z)
end
axis([-300 10 -5 5 -5 5])
hold off

%{
f3=figure;
grid on, xlabel t;    
hold on;
for k=1:1
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    x = matrix(:,2);    
    px = matrix(:,5);
    plot(t,x,'b',t,px,'r');
    legend('x','px');
end
hold off;
%}

nexttile;
grid on, xlabel t, ylabel 'L_x';    
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    x = matrix(:,2);
    y = matrix(:,3);
    z = matrix(:,4);
    px = matrix(:,5);
    py = matrix(:,6);
    pz = matrix(:,7);
    lx = y.*pz-z.*py;
    ly = z.*px-x.*pz;
    lz = x.*py-y.*px;
    %plot(t,lx,t,ly,t,lz);
    %legend('Lx','Ly','Lz');
    plot(t,lx);
    
end
hold off;


saveas(f,'plot.jpg');









function [matrix] = ReadTxt(path, head)
    file = fopen(path, 'r');
    lts=0;
    if (head) lts=1; end
    matrix = dlmread(path,'',lts,0);
    %nPoints = size(matrix,1);
    %nCols = size(matrix,2);
    fclose(file);
end