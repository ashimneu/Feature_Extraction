function [fig] = plot_clusters(fignum,pc_xyz,cluster_ids,pc_idx)
    if nargin == 4
        input_idx = pc_idx;
    else
        input_idx = [];
    end

    [uniq_id,~,~] = unique(cluster_ids);
    
    % mrkr_colors = ['b','r','k','g','c','m','b']'; % scatter plot marker colors
    mrkr_colors = {[0,0,0],[.3,.75,.93],[.64,.08,.18],[0,.45,0.74],[.5,.18,.56]};
    mrkr_sz = 25; % scatter plot marker size
    total_colors = numel(mrkr_colors);

    num_clust = numel(uniq_id);    % total clusters (including outliers)
    lgnd_list = cell(num_clust,1); % allocate variable for legend entries
    
    % Plot all clusters iteratively
    f = figure(fignum); 
    clf; hold on; grid on

    for j = 1:num_clust
        cluster_id = uniq_id(j);
        idx_id = find(cluster_ids == cluster_id); % index of id-th cluster points
        clust_xyz = pc_xyz(idx_id,:); % jth cluster points
        num_pnt = numel(idx_id); % number of points in jth cluster
        if cluster_id == -1
            % assign green color for outlier points only
            color_rgb = [0,1,0]; % green
        else
            color_num = mod(j,total_colors)+1;  % color number
            color_rgb = mrkr_colors{color_num}; % selected color
        end
        sc = scatter3(clust_xyz(:,1),clust_xyz(:,2),clust_xyz(:,3),mrkr_sz,'.','MarkerEdgeColor',color_rgb);
        if ~isempty(input_idx)
            tip_clust_i_idx = input_idx(idx_id); 
        else 
            tip_clust_i_idx = idx_id;  
        end
        sc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',tip_clust_i_idx); 
        tip_clust_i_name = cluster_id.*ones(num_pnt,1); % values for 'cluster' data tip
        sc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cluster',tip_clust_i_name);    
        lgnd_list{j} = strcat(num2str(cluster_id));
    end    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend(lgnd_list,'Location','eastoutside')
    title('Clusters')
    axis equal
    view([15 30])
    
    % return figure handle if output argument is provided
    if nargout == 1        
        fig = f;
    end

end

