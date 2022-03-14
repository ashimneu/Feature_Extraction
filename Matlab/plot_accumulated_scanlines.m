% Please run this script after running process_clustered_scanline_v2.m
% Required variable from process_clustered_scanline_v2.m: scanlist

clc;
addpath(genpath('./utilities'));

% Loading data to matlab workspace takes a few seconds. 
% So, run this line only once.
if ~exist('scanlist','var')
    load('./data/scanlist.mat'); 
end

% parameters: Misc.
p.lasernum = 1:8; % laser num [1 through 8] for feature processing
% parameters: Figure
p.mrk_size = 3;  % default marker size for uncategorized points
p.pnt_size = 20; % Plotting size of feature points
p.s = rng; % random number generator

for scanid = 1:50 % 101 %1:30:length(list) 
    clearaxis(2999); % figure of accumulated scanlines of previous scan is cleared

    fprintf('\n scanid: %3.0f \n',scanid); % print current Lidar scan number
   
    for lzid = p.lasernum % laser_ID
        pc_xyz = scanlist(scanid).laser(lzid).scn_ln.Location;
        save('./utilities/pc_xyz.mat','pc_xyz'); % for distance compute using 'disti.m'
        pc_xyz(:,3) = -pc_xyz(:,3); % invert z-vals for better visualization
        range  = scanlist(scanid).laser(lzid).range;
        theta  = scanlist(scanid).laser(lzid).theta;
        f = scanlist(scanid).laser(lzid).f;
        cluster_ids = scanlist(scanid).laser(lzid).cluster_ids;
        idx_c = scanlist(scanid).laser(lzid).idx_c; % indices of corner points
        idx_nc = scanlist(scanid).laser(lzid).idx_nc; % indices of non-corner points
        pc_c = pc_xyz(idx_c,:); % corner points 
        pc_nc = pc_xyz(idx_nc,:); % non-corner points        
          
        % Plot features for points from scanline of laser number = lzid
        fp = p;
        fp.eb_feats = false; 
        fp.eb_clusteronly = false;
        plot_geometric_features_v3(fp,300,f,lzid);
        
        % plot laser scanline & highlight corner, non-corner points
        plot_scanline(2999,pc_xyz,cluster_ids,idx_c,idx_nc,pc_c,pc_nc,f,lzid);
        title({['Scan #',num2str(scanid),', accumulated scanlines']})

        % hist_point_separation(pc_xyz,10000); % histogram for point separation visualization
        % hist_rangetheta_separation(range,theta,10000); % histogram for range & theta visualization
        % plot_clusters(999,pc_xyz,cluster_ids); % Plot clusters generated on 1 laser scanline
        % title({['scanline, laser',num2str(lzid)]})
    end
end