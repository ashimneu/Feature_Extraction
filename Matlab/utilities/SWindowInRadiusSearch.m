function [pc_sw_out,idx_pnt_ref_sw,idx_out] = SWindowInRadiusSearch(radius,pc,idx_ref_pc)
    % radius - radius of sphere about which neighboring points are selected.
    % pc     - [num x 3] set of points in a scanline
    % idx_ref_pc - index of reference point about in the pc array for which
    %              neighboring points are selected
    % pc_sw  - points selected as neighbors using the sliding window.
    % idx_pnt_ref_sw - index of reference point in the sliding window.
    % index_out - selected neighbors index in pc array.
    pnt_ref = pc(idx_ref_pc,:);
    pntcloud = pointCloud(pc);
    [idx_nei,~] = findNeighborsInRadius(pntcloud,pnt_ref,radius);
    [sorted_idx_nei,sorting_idx] = sort_index(idx_nei);
    pc_sw_out = pntcloud.Location(sorted_idx_nei,:);
    % pc_sw_out = pc_sw(sorting_idx,:);
    
    % index of ref in the array of the neighborhood
    idx_pnt_ref_sw = find(sorted_idx_nei == idx_ref_pc); 

    if nargout == 3
        idx_out = sorted_idx_nei;
    end

    DEBUG = false;
    if DEBUG
        figure(99); clf; hold on; grid on
        pc_sw3 = pntcloud.Location(idx_nei,:);
        % pc_sw3 = pc_sw_out; 
        % pc_sw3 = pc_sw;
        plt = plot3(pc_sw3(:,1),pc_sw3(:,2),pc_sw3(:,3),'.');
        indices1 = 1:size(pc_sw3,1);
        plt.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',indices1);
    end
end

