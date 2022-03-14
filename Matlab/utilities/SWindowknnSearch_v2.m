function [pc_sw,idx_pnt_ref_sw,index_out] = SWindowknnSearch_v2(pc,idx_ref,hf_sz)
    % Select points of pc using a sliding window. 
    % knn neigbhors are selected for a window of size k, centered it ith
    % point of the pc. k should be an odd integer. half of the points
    % should be preceeding point and halft should be succedding point from
    % the ith point.
    % pc   - point cloud
    % idx_ref - index of reference point about in the pc array for which
    %           neighboring points are selected
    % pc_sw - points selected as neighbors using the sliding window.
    % idx_pnt_ref_sw - index of reference point in the sliding window.
    %                  It's the central point in SW except when the leading
    %                  & trailing ends are being considered.
    % index_out - index of selected neighbors in pc array.

    
    EXCLUDE_SW_FOR_ENDS = false;

    num_pnts = size(pc,1);
    sz_window = 2*hf_sz + 1;
    % The following if-statement considers the case arising at the begining or end of a scanline.
    start_idx = idx_ref - hf_sz; % index of initial point in sliding window
    stop_idx  = idx_ref + hf_sz; % index of last point in sliding window
    pc_sw = [];
    if (start_idx <= 0) && ~EXCLUDE_SW_FOR_ENDS % check if start index is at leading end of the scanline
        start_idx = 1;
        stop_idx  = sz_window;
        pc_sw = pc(start_idx:stop_idx,:);
        idx_pnt_ref_sw = idx_ref;
    elseif (stop_idx > num_pnts) && ~EXCLUDE_SW_FOR_ENDS % check if stop index is at trailing end of the scanline
        start_idx = num_pnts - sz_window + 1;
        stop_idx  = num_pnts;
        pc_sw = pc(start_idx:stop_idx,:); 
        idx_pnt_ref_sw = (hf_sz + 1) + ((idx_ref + hf_sz) - num_pnts);
    elseif  (start_idx > 0) && (stop_idx <= num_pnts)
        start_idx = idx_ref - hf_sz;
        stop_idx  = idx_ref + hf_sz;
        pc_sw = pc(start_idx:stop_idx,:); 
        idx_pnt_ref_sw = hf_sz + 1; % mid point of sw
    end
    
    if nargout == 3   
        index_out = start_idx:stop_idx;
    end

end

