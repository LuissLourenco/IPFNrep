close all 


pasta = '../outputs/';
Files = dir(pasta);
Names = {Files.name};
Names = Names(3:end);
Names = string(Names);


Names = cell(1,1);
Names{1} = sprintf('meme/Out00000.txt');
Names = string(Names);


%{
Names = cell(1,2);
Names{1} = sprintf('Data01_a0_20_p0_10/Out00000.txt');
Names{2} = sprintf('Data0-1_a0_20_p0_10/Out00000.txt');

%Names{1} = sprintf('teste3/Out00000.txt');
%Names{2} = sprintf('Data01_a0_10_p0_25/Out00010.txt');

%Names{1} = sprintf('Teste_pl02/Out00000.txt');
Names = string(Names);
%}


n = size(Names,2);
res = cell(n,1);
n_points =1;

for k=1:n
    [matrix,n_points] = ReadTxt(strcat(pasta,Names(k)),true);
    n_points = size(matrix,1);
    %matrix = matrix(1:n_points,:);
    res{k} = mat2cell(matrix, size(matrix,1), size(matrix,2));
end

azul = [0, 0.616, 0.878];
laranja = [0.922, 0.647, 0.02];
vermelho = [0.784, 0.145, 0.014];
cores = [vermelho; laranja];

cores1 = [vermelho];
for i=1:100 
    cores1 = [cores1; vermelho*(100-i)/100+laranja*i/100];
end
for i=1:100 
    cores1 = [cores1; laranja*(100-i)/100+azul*i/100];
end
for i=1:100 
    cores1 = [cores1; azul*(100-i)/100+vermelho*i/100];
end

cores = [cores1; cores1]

leg = cell(1,2);
leg{1} = sprintf('off');
leg{2} = sprintf('on');
leg = string(leg);


f=figure;
view(axes(),3);
set(gca,'FontSize',12)
set(gcf, 'Position',  [100, 100, 1000, 1000])
%grid on, xlabel 'x (\lambda)', ylabel 'y (\lambda)', zlabel 'z (\lambda)';
grid on;
%title trajectory;
title({''});
set(gca,'visible','off')
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    n_points = n_points/1.9;
    start = 1000;
    t = [matrix(start:n_points,1); NaN];
    x = [matrix(start:n_points,2); NaN];
    y = [matrix(start:n_points,3); NaN];
    z = [matrix(start:n_points,4); NaN];
    
    cla
    meme = patch(x,y,z,x,'EdgeColor','interp','FaceColor','none');
    
    %plot3(x,y,z,'DisplayName',leg(k),'Color',cores(k+2,:));
end

colormap (cores)
%colormap(f);
%legend()
axis([-12000 1000 -20 7 -20 7])

%{
f = figure;
set(gca,'FontSize',12);
grid on, xlabel 'x (\lambda)', ylabel 'y (\lambda)';
hold on;
title({'Linear Polarization'}, {'counter-propagation'});
for k=1:n
    matrix = cell2mat(res{k});
    x = matrix(:,2);
    y = matrix(:,3);       
    plot(x,y,'Color',cores(k+1,:));
end
axis([-18 -12 -1.5 1.5])
%}
nc = 2;
nl = 3;

%{
f2=figure;
set(gcf, 'Position',  [100, 100, 1800, 800])
subplot(nc,nl,1);
set(gca,'FontSize',12)
grid on, xlabel 't (T)', ylabel 'x (\lambda)';
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    x = matrix(:,2);
    plot(t,x,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))),...
        'Color',cores(k,:));
end
ylim([-1000 1000])
%legend()

subplot(nc,nl,2);
set(gca,'FontSize',12)
%set(gca, 'YScale', 'log');
grid on, xlabel 't (T)', ylabel 'y (\lambda)';
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    y = matrix(:,3);
    plot(t,y,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))),...
        'Color',cores(k,:));
        
end
ylim([-10 10])
%legend()

subplot(nc,nl,3);
set(gca,'FontSize',12)
grid on, xlabel 't (T)', ylabel 'z (\lambda)';
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    z = matrix(:,4);
    plot(t,z,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))),...
        'Color',cores(k,:));
end
ylim([-10 10])
%legend()


subplot(nc,nl,4);
grid on, xlabel t, ylabel px;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    px = matrix(:,5);
    plot(t,px,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))),...
        'Color',cores(k,:));
end
%legend()

subplot(nc,nl,5);
grid on, xlabel t, ylabel py;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    y = matrix(:,3);
    z = matrix(:,4);
    py = matrix(:,6);
    pz = matrix(:,7);
    phi = atan2(z,y);
    pr = py.*cos(phi) + pz.*sin(phi);
    plot(t,py,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))),...
        'Color',cores(k,:));
end
%legend()

subplot(nc,nl,6);
grid on, xlabel t, ylabel pz;
hold on;
for k=1:n
    matrix = cell2mat(res{k});
    t = matrix(:,1);
    y = matrix(:,3);
    z = matrix(:,4);
    py = matrix(:,6);
    pz = matrix(:,7);
    phi = atan2(z,y);
    pphi = -py.*sin(phi) + pz.*cos(phi);
    plot(t,pz,'DisplayName',strcat('\phi_0 = ',num2str(phi(1))),...
        'Color',cores(k,:));
end
%legend()
%}

 






function [matrix, nPoints] = ReadTxt(path, head)
    file = fopen(path, 'r');
    lts=0;
    if (head) lts=1; end
    matrix = dlmread(path,'',lts,0);
    nPoints = size(matrix,1);
    %nCols = size(matrix,2);
    fclose(file);
end