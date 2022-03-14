% This is the final version submitted at the end of the internship.
% Author: Ashim Neupane, email: aneup001@ucr.edu

clear; clc;
addpath(genpath('./utilities'));
[path,scanlist,datasetprefix] = loadDatasetInfo(3);
addpath(genpath(path));
scanlist = sortbyfilenumber(scanlist,datasetprefix);  % sort struct entries in list by number in filename
laserID = load_laserID(path,datasetprefix);

% parameters: Neighbors Selection
p.sz_window = 21; % window size 
p.sw_def = 2; % neighbors selection: 1, 2 & 3 - knnSearch v1, v2, sphere.
p.radius = 1; % [meters] radius of sphere for neighboring point selection
% parameters: Clustering Thresholds 
p.range_threshold = 0.5; % [meters]
p.theta_threshold = 0.003; % [degrees]
p.minpnts = 10;
% parameters: Feature Thresholds
p.threshold_Sm = 0.005; % 0.04
p.threshold_L = 0.85; % for extraction of corner use 0.7, for planes use 0.9
p.threshold_dp = 0.6;   % dotproduct
% parameters: Misc.
p.lasernum = 1:8; % laser num [1 through 8] for feature processing
% parameters: Figure
p.mrk_size = 3;  % default marker size for uncategorized points
p.pnt_size = 20; % Plotting size of feature points
p.s = rng; % random number generator

% declare variables to store pointclouds.
pc_c = []; pc_l = []; pc_all = [];

for scanid = 1:900 % 101 %1:30:length(list) 
    % clearfig(2999); % figure of accumulated scanlines of previous scan is cleared

    fprintf('\n scanid: %3.0f \n',scanid); % print current Lidar scan number

    % load dataset of current scan
    scanlist(scanid).pc = pcread(scanlist(scanid).folder + "\" + scanlist(scanid).name);
    pc_all = [pc_all; scanlist(scanid).pc.Location];

    % Extract & group together points of each laser
    [scanlist(scanid).laser] = scanline_pc_extract(scanlist(scanid).pc,laserID(scanid)); 
    
    for lzid = p.lasernum % laser_ID
        pc_xyz = scanlist(scanid).laser(lzid).scn_ln.Location;
        save('./utilities/pc_xyz.mat','pc_xyz'); % for distance compute using 'disti.m'
        pc_xyz(:,3) = -pc_xyz(:,3); % invert z-vals for better visualization
        range  = scanlist(scanid).laser(lzid).range;
        theta  = scanlist(scanid).laser(lzid).theta;
        % pc_index = 1:size(pc_xyz,1); % indices of points in pc_xyz array

        % clustering of one laser scanline based on differences in range & theta values       
        cluster_ids = clustering_using_rangethetav2(range,theta,p.range_threshold,p.theta_threshold,p.minpnts);

        % hist_point_separation(pc_xyz,10000); % histogram for point separation visualization
        % hist_rangetheta_separation(range,theta,10000); % histogram for range & theta visualization
        % plot_clusters(999,pc_xyz,cluster_ids); % Plot clusters generated on 1 laser scanline
        % title({['scanline, laser',num2str(lzid)]})

        idx_endpoint = size(pc_xyz,1); % index of last data in the scanline
        num_pnts     = size(pc_xyz,1); % number of points in the scanline
        
        [uniq_cid,~,~] = unique(cluster_ids); % count of unique cluster ids

        % Declare variables to store computed feature values. 
        f = create_feature_struct(num_pnts);
        avg_density = zeros(max(uniq_cid),1);
        avg_dist  = zeros(max(uniq_cid),1); % [meter] average point separation in cluster id
        pntdist   = nan(num_pnts,1); % distance between points whose indices are i & i+1
        ad_radius = zeros(max(uniq_cid),1); % [meter] adaptive radius for cluster id
        
        for cluster_id = 1:max(uniq_cid)
            idx_clust_id = find(cluster_ids == cluster_id); % index of id-th cluster points in the scanline (i.e. pc_xyz) 
            pnt_clust = pc_xyz(idx_clust_id,:); % points in id-th cluster
            
            % check if cluster contains leading & trailing ends of the scanline.
            % if yes, to avoid index dropping from end to 1 during feature
            % computation, for points 
            % whose index = 1, new index = last index + 1, 
            % whose index = 2, new index = last index + 2, ....
            if  contains_ends(idx_clust_id,idx_endpoint)
                idx_clust_id2 = shift_index(idx_clust_id);
                [idx_clust_id2,sorting_idx] = sort(idx_clust_id2,'ascend');
                pnt_clust2   = pnt_clust(sorting_idx,:);
                idx_clust_id = idx_clust_id(sorting_idx);
            else
                idx_clust_id2 = idx_clust_id;
                pnt_clust2    = pnt_clust;
            end
        
            num_pnt_id = size(pnt_clust,1); % number of points in idth cluster
            % sz_windows = 5:2:odd_int(num_pnt); list of window sizes to evalute
            
            % compute density over the cluster points
            count_pnt = 3; % number of points over which density is computed
            [density(idx_clust_id),avg_dist(cluster_id),pntdist(idx_clust_id)] = compute_density(pnt_clust2,count_pnt); % [# points/meter]
            avg_density(cluster_id) = mean(density(idx_clust_id),'omitnan'); % average density in the cluster
            pntsep_prctile = prctile(pntdist(idx_clust_id),95); % pnt separation at selected percentile
        
            % Compute adaptive radius for sphere based neighbor selection method
            neigbhorhood_pnt_count = 10;
            ad_radius(cluster_id)  = neigbhorhood_pnt_count*pntsep_prctile;
            p.radius = ad_radius(cluster_id);
        
            p.idxa = idx_clust_id; % DELETE ME AFTER DEBUGGING!
            % compute geometric features on idth cluster
            f_cluster_id = compute_feature_using_sw_v4(p,pnt_clust2,p.sz_window);  
        
            % % Plot features for current cluster
            % fp = p; 
            % fp.eb_feats = false; % enable/disable this plot
            % fp.eb_clusteronly = true;
            % fp.sz_window = p.sz_window; fp.clust_id = cluster_id;
            % plot_geometric_features(fp,301,features_cluster_id,lzid);
            
            % Save computed features
            f = fill_feature_struct(f,f_cluster_id,idx_clust_id);
        
            % h = hist_point_separation(pnt_clust,10000); % histogram for point separation visualization
            % title({['scanline, laser',num2str(lzid),', cluster ',num2str(cluster_id)]})
        end
    
        % Plot features for all clusters (except outliers)
        fp = p; 
        fp.eb_feats = false; 
        fp.eb_clusteronly = false;
        plot_geometric_features_v3(fp,300,f,lzid);
        
        %% Accumulate & plot scanline of multiple lasers
        % threshold check on feature: Linearity/Smoothness/Dotproduct
        % [idx_g,idx_s,pc_g,pc_s] = find_index(f.Sm,p.threshold_Sm,pc_xyz);
        % [idx_s,idx_g,pc_s,pc_g] = find_index(f.L,p.threshold_L,pc_xyz);
        [idx_s,idx_g,pc_s,pc_g] = find_index(f.dotprod,p.threshold_dp,pc_xyz);

        % save cluster_ids, features & indices of corner & non-corner points
        scanlist(scanid).laser(lzid).f = f;
        scanlist(scanid).laser(lzid).cluster_ids = cluster_ids;
        scanlist(scanid).laser(lzid).idx_c  = idx_s;
        scanlist(scanid).laser(lzid).idx_nc = idx_g;
        
        % % plot laser scanline & highlight corner, non-corner points
        % plot_scanline(2999,pc_xyz,cluster_ids,idx_s,idx_g,pc_s,pc_g,f,lzid);
        % title({['Scan #',num2str(scanid),', accumulated scanlines']})
    end
end

% save('./data/scanlist.mat','scanlist')
% save('./data/scanlist.mat','scanlist','-v7.3')




