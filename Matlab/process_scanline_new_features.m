% This script is a newer variant of Zeyi Jiang's process_scanline.m.
% Thanks to Zeyi Jiang for providing initial iteration of process_scanline.m 
% and its supporting scripts.
% process_scanline_new_feautres.m was created to avoid conflicts during 
% commit while other newer geometric features are being test on lidar 
% scanlines for Microstar's Container Project.
% Author: Ashim Neupane, email: aneup001@ucr.edu

clear; clc;
addpath('./Test_Features'); % for compute_moment2.m
addpath(genpath('./utilities'));
[path,list,datasetprefix] = loadDatasetInfo(3);
addpath(genpath(path));
list = sortbyfilenumber(list,datasetprefix);  % sort struct entries in list by number in filename
laserID = load_laserID(path,datasetprefix);

% SW parameters
p.sz_window = 15; % window size (use 151 for straight line detection & 351 for corner detection; works well on scan 600, laser 6)
p.sw_def = 3; % 1,2,3 are knnSearch v1,v2 and distance based window respectively 
p.radius = 5; % 1.75; % [meters] sphere radius to define sliding window neighborhood 
% threshold parameters
p.range_threshold = 0.5;
p.theta_threshold = 0.01;
p.threshold_Sm = 0.04; % 0.04
p.threshold_L = 0.85; % for extraction of corner use 0.7, for planes use 0.9
% misc. parameters
p.lasernum = 3; % [1 through 8] specify laser num for feature pc processing
p.dcsnfeat_cor = 1; % decision feature(s):, 1 - L, 2 - Sm, 3 - L & Sm
% figure parameters
p.eb_moment = 1; % 1-enable, 0-disable, moment figure
p.eb_sm = 1;     % 1-enable, 0-disable, smoothness figure
p.eb_feats = 1;  % 1-enable, 0-disable, the features: linearity, Surf. Var., Omnivar, Eigen Entropy
p.mrk_size = 3;  % default marker size for uncategorized points
p.pnt_size = 20; % Plotting size of feature points
p.s = rng;

% declare variables to store pointclouds.
pc_c = []; pc_l = []; pc_all = [];

for scanid = 1000 %1:30:length(list)    
    fprintf('\n'); fprintf('i %2.0f \n',scanid);

    % load dataset of current scan
    list(scanid).pc = pcread(list(scanid).folder + "\" + list(scanid).name);
    pc_all = [pc_all; list(scanid).pc.Location];

    % Extract & group together points of each laser
    [list(scanid).laser] = scanline_pc_extract(list(scanid).pc,laserID(scanid)); 
    
    % f1 = figure(1000); clf; hold on; grid on    

    for lzid = p.lasernum % laser_ID
        pc_xyz   = list(scanid).laser(lzid).scn_ln.Location;
        pc_xyz(:,3) = -pc_xyz(:,3); % invert z-vals for better visualization
        range    = list(scanid).laser(lzid).range;
        theta    = list(scanid).laser(lzid).theta;
        pc_index = 1:size(pc_xyz,1); % indices of points in pc_xyz array

        % % histogram for point separation visualization
        % hist_point_separation(pc_xyz,10000);

        % clustering of one laser scanline based on differences in range & theta values       
        cluster_ids = clustering_using_rangetheta(p,range,theta);
        
        % % Draw figure for the cluster visualization            
        % plot_clusters(999,pc_xyz,cluster_ids); % fignum:999
        % title({['pc clustered using range & theta, laser',num2str(lzid)]})

        lazer = struct;
        if (~isempty(list(scanid).laser(lzid).pnt_idx))
           [lazer.cor_pc, lazer.ln_pc, lazer.cor_pc_index,...
               lazer.ln_pc_index, lazer.otherfeatures] = ...
               scanline_feature_extract_v2(p,pc_xyz);
        end
        
        % merge feature struct of multiple lazers for different scanids .
        % See function './utilities/mergestruct.m' for more details.
        fn = fieldnames(lazer);
        for i = 1:length(fn)
            list(scanid).laser(lzid).(fn{i}) = lazer.(fn{i});
        end  
        
        % retrieve computed features for plotting
        cor_pc = list(scanid).laser(lzid).cor_pc(:,1:3); % pc of points detected in corner region
        cor_index = list(scanid).laser(lzid).cor_pc_index;
        ln_pc  = list(scanid).laser(lzid).ln_pc(:,1:3);  % pc of points detected as line points
        ln_index  = list(scanid).laser(lzid).ln_pc_index;
        otherfeatures = list(scanid).laser(lzid).otherfeatures;

        % Plot scanlines
        figure(1000); clf; hold on; grid on    

        % Plot original scanline        
        s1 = scatter3(pc_xyz(:,1),pc_xyz(:,2),pc_xyz(:,3),p.mrk_size,'b');              
        s1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',pc_index);

        % plot points extracted as being corner points
        if(~isempty(list(scanid).laser(lzid).cor_pc))            
            s2 = scatter3(cor_pc(:,1),cor_pc(:,2),cor_pc(:,3),p.pnt_size,'rp');
            s2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',cor_index);
            pc_c = [pc_c;cor_pc];
        end 
        % plot points extracted as being line points
        if(~isempty(list(scanid).laser(lzid).ln_pc))            
            s3 = scatter3(ln_pc(:,1),ln_pc(:,2),ln_pc(:,3),p.pnt_size,'c');
            s3.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',ln_index);
            pc_l = [pc_l;ln_pc];
        end

        % hold off
        title({['scanline, laser',num2str(lzid)]})
        xlabel('x, m')
        ylabel('y, m')
        zlabel('z, m') 
        legend('scnln','cr','non cr')
        % axis equal
        view([0 75])
        xlim([-100 100])
        ylim([-60 40])
%         zlim([-25 0])

        % Do smoothing & Plot current laser scanline
        % sm_window = 10; % smoothing window
        % sm_pc_xyz = doSmoothing_Plot(pc_xyz,sm_window); % smoothens data & does scatter3 plot
        % [dx,dy,dz] = compute_derivative(pc_xyz,'drawfig');
       
        % Plot x,y,z components of points for each laser scanline
        % plot_components(p,1,pc_xyz,lzid); % fig number - 1

        % Plot Geometric Features
        fp = p; fp.eb_clusteronly = false; fp.eb_feats = true;
        plot_geometric_features(fp,100,otherfeatures,lzid);

        % Plot range & theta values of current scanline points
        % plot_range_theta(p,101,range,theta,lzid); % fig number - 101               
    end
    % pc_c = [pc_c;list(i).cor_pc(:,1:3)];
    % pc_l = [pc_l;list(i).ln_pc(:,1:3)];
end

% figure(2); clf %  A canvas for drawing scans at one place
% title('Accumulated scans')
% hold on
% pcshow(pc_all,'c')
% % pcshow(pc_c,'r','MarkerSize',feature_pnt_size)
% % pcshow(pc_l,'y','MarkerSize',feature_pnt_size)
% hold off
% xlabel('x, m')
% ylabel('y, m')
% zlabel('z, m')
% view([0 75])
% xlim([-100 100])
% ylim([-60 40])
% zlim([-25 0])