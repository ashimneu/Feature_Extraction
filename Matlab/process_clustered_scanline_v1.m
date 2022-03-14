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
p.threshold_Sm = 0.04; % 0.04
p.threshold_L = 0.85; % for extraction of corner use 0.7, for planes use 0.9
% parameters: Misc.
p.lasernum = 4; % [1 through 8] specify laser num for feature pc processing
% parameters: Figure
p.mrk_size = 3;  % default marker size for uncategorized points
p.pnt_size = 20; % Plotting size of feature points
p.s = rng;

% declare variables to store pointclouds.
pc_c = []; pc_l = []; pc_all = [];

for scanid = 100 %1:30:length(list)    
    fprintf('\n'); fprintf('i %2.0f \n',scanid);

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

        % % histogram for point separation visualization
        % hist_point_separation(pc_xyz,10000);

        % % histogram for range & theta visualization
        % hist_rangetheta_separation(range,theta,10000);        
    end
end

% % Plot Clusters generate on 1 laser scanline 
% plot_clusters(999,pc_xyz,cluster_ids); % fignum:999
% title({['scanline, laser',num2str(lzid)]})

idx_endpoint = size(pc_xyz,1); % index of last data in the scanline
num_pnts     = size(pc_xyz,1); % number of points in the scanline

[uniq_id,ia,ic] = unique(cluster_ids); % count of all unique clusters

% Declare variables to store computed feature values.
[L,P,S,O,A,Esum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dotprod,len_pr,len_sc,density,angle] = deal(nan(num_pnts,1));
avg_density = zeros(max(uniq_id),1);
avg_dist  = zeros(max(uniq_id),1); % [meter] average point separation in cluster id
pntdist   = nan(num_pnts,1); % distance between points whose indices are i & i+1
ad_radius = zeros(max(uniq_id),1); % [meter] adaptive radius for cluster id

for cluster_id = 1:max(uniq_id)
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
    features_cluster_id = compute_feature_using_sw_v3(p,pnt_clust2,p.sz_window);  

    % Plot features for current cluster
    fp = p; 
    fp.eb_feats = false; % enable/disable this plot
    fp.eb_clusteronly = true;
    fp.sz_window = p.sz_window; fp.clust_id = cluster_id;
    plot_geometric_features(fp,300,features_cluster_id,lzid);
    
    % Save computed features
    [L(idx_clust_id),P(idx_clust_id),S(idx_clust_id),O(idx_clust_id),A(idx_clust_id),Esum(idx_clust_id),SV(idx_clust_id),...
    Entropy(idx_clust_id),Sm(idx_clust_id),FOFA(idx_clust_id),FOSA(idx_clust_id),SOFA(idx_clust_id),...
    SOSA(idx_clust_id),dotprod(idx_clust_id),angle(idx_clust_id)] = deal(features_cluster_id{:});
    
    % % histogram for point separation visualization
    % h = hist_point_separation(pnt_clust,10000);
    % title({['scanline, laser',num2str(lzid),', cluster ',num2str(cluster_id)]})
end

features_allcluster = {L,P,S,O,A,Esum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dotprod,angle};

% Plot features for all clusters (except outliers)
fp = p; 
fp.eb_feats = true; 
fp.eb_clusteronly = false;
plot_geometric_features_v2(fp,300,features_allcluster,lzid);

% % Plot density for all clusters (except outliers)
% figure(399); clf; hold on; grid on
% plot(density,'.')
% xlabel('index'); ylabel('density, count/meter')

%% Accumulate & plot scanline of multiple lasers

% parameters: Feature Thresholds
p.threshold_Sm = 0.005; % smoothness
p.threshold_L  = 0.85;  % linearity 
p.threshold_dp = 0.6;   % dotproduct

% enable only 1 of following three threshold check on Linearity, Smoothness & dotproduct.
% [idx_g,idx_s,pc_g,pc_s] = find_index(Sm,p.threshold_Sm,pc_xyz);
% [idx_s,idx_g,pc_s,pc_g] = find_index(L,p.threshold_L,pc_xyz);
[idx_s,idx_g,pc_s,pc_g] = find_index(dotprod,p.threshold_dp,pc_xyz);

% features struct
features.angle   = angle;
features.dotprod = dotprod;
plot_scanline(2999,pc_xyz,cluster_ids,idx_s,idx_g,pc_s,pc_g,features);


% % Parameters: figure
% sz_mrk = 75; % marker size
% 
% figure(2999); clf; hold on; grid on;
% % all points
% plt_all = scatter3(pc_xyz(:,1),pc_xyz(:,2),pc_xyz(:,3),sz_mrk,'.','MarkerEdgeColor',[0,1,0]); % corner
% plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',pc_index);
% plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cluster',cluster_ids);
% plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dotprod',dotprod);
% plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('angle',angle);
% % corner points
% plt_c = scatter3(pc_s(:,1),pc_s(:,2),pc_s(:,3),sz_mrk,'.','MarkerEdgeColor',[1,0,0]); % corner
% plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idx_s);
% plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cluster',cluster_ids(idx_s));
% plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dotprod',dotprod(idx_s));
% plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('angle',angle(idx_s));
% % non corner points
% plt_nc = scatter3(pc_g(:,1),pc_g(:,2),pc_g(:,3),sz_mrk,'.','MarkerEdgeColor',[.3,.75,.93]); % non-corner
% plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idx_g);
% plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cluster',cluster_ids(idx_g));
% plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dotprod',dotprod(idx_g));
% plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('angle',angle(idx_g));
% xlabel('x'); ylabel('y'); zlabel('z')
% legend('outliers','cor.','non cor.')
% view([20 30])
% axis equal

% cor_pc = pc_s;
% ln_pc  = pc_g;
% otherfeatures = features_cluster_all;
