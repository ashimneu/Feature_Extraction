clc
clear
% path = 'ess_scans/';
% path = 'data/';

list = dir('data/*.pcd');   % path of the source folder of scans
flag_DEBUG = 1;             % Flag to turn on/off some debugging code
default_size = 6;           % Plotting size of uncategorized points
feature_pnt_size = 20;      % Plotting size of feature poitns
laserID_A = laserID_logParser('data/laserID_log_A.txt');
laserID_P = laserID_logParser('data/laserID_log_P.txt');


pc_c = [];
pc_l = [];
pc_all = [];

for i= 1:length(list)
    i
    list(i).pc = pcread(list(i).folder+"\"+list(i).name);
    pc_all = [pc_all; list(i).pc.Location];
    %     [list(i).cor_pc, list(i).ln_pc]= scanline_feature_extract_m8(list(i).pc, laserID(i).ID); % Processing each scan to extract feature points. (Not very satisfying yet)
    list(i).feature = scanline_feature_extract_m8(list(i).pc, laserID_P(i).ID); % Processing each scan to extract feature points. (Not very satisfying yet)
    
    figure(1)
    hold on
    scatter3(list(i).pc.Location(:,1),list(i).pc.Location(:,2),list(i).pc.Location(:,3),default_size,'b')
    for j = 1:8
        if(~isempty(list(i).feature(j).cor_pc))
            scatter3(list(i).feature(j).cor_pc(:,1),list(i).feature(j).cor_pc(:,2),list(i).feature(j).cor_pc(:,3),feature_pnt_size,'r')
            pc_c = [pc_c;list(i).feature(j).cor_pc(:,1:3)];
        end
        if(~isempty(list(i).feature(j).cor_pc))
            scatter3(list(i).feature(j).ln_pc(:,1),list(i).feature(j).ln_pc(:,2),list(i).feature(j).ln_pc(:,3),feature_pnt_size,'y')
            pc_l = [pc_l;list(i).feature(j).ln_pc(:,1:3)];
        end
        
    end
    hold off
    title('scan display')
    xlabel('x, m')
    ylabel('y, m')
    zlabel('z, m')
    axis equal
end

figure(2)                   %  A canvas for drawing scans at one place
title('Accumulated scans')
hold on
pcshow(pc_all,'b')
pcshow(pc_c,'r','MarkerSize',feature_pnt_size)
pcshow(pc_l,'y','MarkerSize',feature_pnt_size)
hold off
xlabel('x, m')
ylabel('y, m')
zlabel('z, m')
axis equal
