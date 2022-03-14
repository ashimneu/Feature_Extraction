function cluster_ids = clustering_using_rangethetav2(range,angle,range_threshold,angle_threshold,minpnts)
    % Here, a scanline point cloud is separated into clusters based on
    % differences in range and angluar values (theta) of the points.
    % If at any point index the range/theta difference between successive 
    % points is greater than a threshold a new cluster is started by 
    % splitting the current cluster at current point. 
    % If any of the splitting indices are similar between the set of 
    % splitting indexes generated in range & theta based clustering, it is
    % made sure that the redundancy is identified & removed.
    % If any cluster has points less than minimum required for a cluster,
    % that group will be denoted as cluster of outliers. 
    % All clusters are assigned id = +ve integerers >= 1, except cluster
    % of outlier points which is denoted by cluster id = -1.
    % Author : Ashim Neupane, email: aneup001@ucr.edu
    % INPUTS : 
    % range - [numx1] range data
    % angle - [numx1] laser horizontal angle data
    % range_threshold - threshold for range difference based cluster splitting
    % theta_threshold - threshold for theta difference based cluster splitting
    % minpnts - minimum number of points required in a cluster
    % OUTPUT : 
    % cluster_ids - [numx1] array of cluster id (including -1 for outliers)
    % cluster_id = -1 for outlier points
    
    num_pnt = numel(range);
    diff_range = range(2:end) - range(1:end-1);
    diff_theta = angle(2:end) - angle(1:end-1);
    
    % do threshold check
    flag_r = abs(diff_range) >= range_threshold;
    flag_t = abs(diff_theta) >= angle_threshold;
    
    index_fl_range = find(flag_r); % index of flags in range values
    index_fl_theta = find(flag_t); % index of flags in theta values
    
    % check & remove any theta based flag if its index matches with range 
    % based flag to avoid cluster breaking at same index during theta
    % difference based clustering.
    [~,i_fl_range,i_fl_theta] = intersect(index_fl_range,index_fl_theta); %#ok<ASGLU> 
    index_fl_theta(i_fl_theta) = [];
    
    [~,dim_max] = max(size(range)); % get dimension number along length of the array
    index_all = cat(dim_max,index_fl_range,index_fl_theta);
    index_all = sort(index_all,'ascend'); 
    num_cluster = numel(index_all)+1;   % +1 because dividing an array using
                                        % n breaking index results in n+1
                                        % children arrays.
    true_cluster_id = 1;
    % Cluster point separation based on breaking indices
    for cluster_id = 1:num_cluster
        if cluster_id == 1
            start_idx = 1;
        else
            start_idx = index_all(cluster_id-1)+1;
        end
        if cluster_id == num_cluster
            end_idx = num_pnt;
        else
            end_idx   = index_all(cluster_id);
        end
        
        % denote a cluster as a group of outlier points if point count is 
        % less than minimum required points.
        num_pnt_clust = end_idx - start_idx + 1; % number of points
        if num_pnt_clust <= minpnts
            cluster_ids(start_idx:end_idx) = -1; % outlier
        else
            cluster_ids(start_idx:end_idx) = true_cluster_id;
            true_cluster_id = true_cluster_id + 1;
        end
    
        % If initial and final clusters are nearby each other and fall
        % within range & theta threshold, fuse them.
        if cluster_id == num_cluster
            range_is_good = abs(range(num_pnt) - range(1)) < range_threshold;
            theta_is_good = abs(angle(num_pnt) - 2*pi - angle(1)) < angle_threshold;
            if range_is_good && theta_is_good
                cluster_ids(start_idx:end_idx) = 1;
            end
        end
    end
end

