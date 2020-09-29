close all 


pasta = 'Outputs/';
Files=dir(pasta);
Names = {Files.name};
Names=Names(3:end);
Names = string(Names);
n = size(Names,2);
res = cell(n,1);

for k=1:n
    [matrix] = ReadTxt(strcat(pasta,Names(k)),false);
    n_points = size(matrix,1);
    %matrix = matrix(1:n_points/2,:);
    res{k} = mat2cell(matrix, size(matrix,1), size(matrix,2));
end



f=figure;
view(axes(),3);
grid on, xlabel x, ylabel y, zlabel z;
hold on;
for k=[10 2]
    matrix = cell2mat(res{k});
    x = matrix(:,2);
    y = matrix(:,3);    
    z = matrix(:,4);
    plot3(x,y,z);
end
%axis([-300 10 -1 1 -1 1])


%{
f2 = figure;
grid on, xlabel y, ylabel z;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    x = matrix(:,2);
    y = matrix(:,3);    
    z = matrix(:,4);
    plot(y,z);
end
hold off;
legend([{'with E_x'},{'without E_x'}]);

f3 = figure;
grid on, xlabel t, ylabel p_r;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    y = matrix(:,3);    
    z = matrix(:,4);
    phi = atan2(z,y);
    py = matrix(:,6);
    pz = matrix(:,7);
    pphi = -sin(phi).*py+cos(phi).*pz;
    pr = py.*cos(phi)+pz.*sin(phi);
    plot(t,pr);
end
hold off;
legend([{'with E_x'},{'without E_x'}]);


f4 = figure;
grid on, xlabel t, ylabel p_x;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    px = matrix(:,5);
    plot(t,px);
end
hold off;
legend([{'with E_x'},{'without E_x'}]);

%}










function [matrix] = ReadTxt(path, head)
    file = fopen(path, 'r');
    lts=0;
    if (head) lts=1; end
    matrix = dlmread(path,'',lts,0);
    %nPoints = size(matrix,1);
    %nCols = size(matrix,2);
    fclose(file);
end