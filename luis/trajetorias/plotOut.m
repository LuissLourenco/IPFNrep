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
    matrix = matrix(1:n_points,:);
    res{k} = mat2cell(matrix, size(matrix,1), size(matrix,2));
end



f=figure;
view(axes(),3);
grid on, xlabel x, ylabel y, zlabel z;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    x = matrix(:,2);
    y = matrix(:,3);    
    z = matrix(:,4);
    phi = atan2(z,y);
    plot3(x,y,z,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))));
end
legend()
%axis([-1800 800 -5 5 -5 5])

f2=figure;
grid on, xlabel t, ylabel y;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    x = matrix(:,2);
    y = matrix(:,3);
    z = matrix(:,4);
    px = matrix(:,5);
    py = matrix(:,6);
    phi = atan2(z,y);
    plot(t,y,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))));
end
legend()
%axis([0 3000 -10 10])

 






function [matrix] = ReadTxt(path, head)
    file = fopen(path, 'r');
    lts=0;
    if (head) lts=1; end
    matrix = dlmread(path,'',lts,0);
    %nPoints = size(matrix,1);
    %nCols = size(matrix,2);
    fclose(file);
end