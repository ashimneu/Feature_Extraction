function f = fill_feature_struct(f,f_cluster_id,idx_clust_id)
    % copies features computed for idth cluster points of a scanline
    % into the larger struct f which accumulates the features for 
    % all the points in the scanline.

    % f - struct of features of points in the scanline
    % f_cluster_id - struct of features of points in idth cluster of the scanline
    % idx_clust_id - index of points in idth cluster in scanline array (i.e. pc_xyz)

    fn_f = fieldnames(f);            % fieldnames of f
    fn_c = fieldnames(f_cluster_id); % fieldnames of f_cluster_id
    cfields = intersect(fn_f,fn_c);  % find common fields
    
    % iterate over the common features in both structs
    for i = 1:numel(cfields)
        cfield = cfields{i};
        f.(cfield)(idx_clust_id') = f_cluster_id.(cfield);
    end
end

