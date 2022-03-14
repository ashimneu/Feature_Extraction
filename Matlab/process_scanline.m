clc
clear

addpath('./Test_Features')
list = dir('data/A/*.pcd');     % path of the source folder of scans
list = sortbyfilenumber(list);  % sort struct entries in list by number in filename
flag_DEBUG   = 1;               % Flag to turn on/off some debugging code
default_size = 6;               % Plotting size of uncategorized points
feature_pnt_size = 20;          % Plotting size of feature poitns


pc_c = [];
pc_l = [];
pc_all = [];

for i= 1:length(list) 
    i
    list(i).pc = pcread(list(i).folder+"\"+list(i).name);
    pc_all = [pc_all; list(i).pc.Location];
    [list(i).cor_pc, list(i).ln_pc]= scanline_feature_extract(list(i).pc); % Processing each scan to extract feature points. (Not very satisfying yet)
    
    % Plot current scanned point cloud
    figure(1)
    scatter3(list(i).pc.Location(:,1),list(i).pc.Location(:,2),list(i).pc.Location(:,3),default_size,'b')
    hold on
    scatter3(list(i).cor_pc(:,1),list(i).cor_pc(:,2),list(i).cor_pc(:,3),feature_pnt_size,'r')
    %     scatter3(list(i).ln_pc(:,1),list(i).ln_pc(:,2),list(i).ln_pc(:,3),feature_pnt_size,'y')
    hold off
    title('scan display')
    xlabel('x, m')
    ylabel('y, m')
    zlabel('z, m')
    axis equal

    pc_c = [pc_c;list(i).cor_pc(:,1:3)];
    pc_l = [pc_l;list(i).ln_pc(:,1:3)];

end

figure(2); clf                   %  A canvas for drawing scans at one place
title('Accumulated scans')
hold on
pcshow(pc_all,'c')
pcshow(pc_c,'r','MarkerSize',feature_pnt_size)
pcshow(pc_l,'y','MarkerSize',feature_pnt_size)
hold off
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')
axis equal
view([0 90])
