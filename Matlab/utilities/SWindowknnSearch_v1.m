function [pc_sw,idx_pnt_ref_sw,index_out] = SWindowknnSearch_v1(pc,i,hf_sz)
    % Select points of pc using a sliding window. 
    % knn neigbhors are selected for a window of size k, centered it ith
    % point of the pc. k should be an odd integer. half of the points
    % should be preceeding point and halft should be succedding point from
    % the ith point.
    % pc   - point cloud
    % i    - index of point about which sliding window for neigbhoring 
    %        points are extracted.
    
    num_pnts = size(pc,1);
    first_pnt_idx = i - hf_sz;  % first point index
    last_pnt_idx = i + hf_sz;   % last point index

    % The following if-statement considers the case at the begining or end of
    % a scan. Basically, it considers a scan as a circular scan. E.g., If 
    % the point chosen is the first point of a scan, then half of the slide 
    % window will be at the begining of the scan, the other half is the 
    % end of the scan.
    if first_pnt_idx <= 0 % at head
        head_idx = num_pnts - abs(first_pnt_idx);
        tail_idx = i+hf_sz;
        head_pc = pc(head_idx:num_pnts,:);
        pc_sw = [head_pc; pc(1:tail_idx,:)];
    elseif last_pnt_idx > num_pnts % at tail
        head_idx = i-hf_sz;
        tail_idx = last_pnt_idx - num_pnts;
        tail_pc = pc(1:tail_idx,:);
        pc_sw = [ pc(head_idx:num_pnts,:);tail_pc];
    else
        head_idx = i-hf_sz;
        tail_idx = i+hf_sz;
        pc_sw =  pc(head_idx:tail_idx,:);
    end
    
    idx_pnt_ref_sw = hf_sz + 1;

    if nargout == 3
        if head_idx > tail_idx
            index_out = [head_idx:num_pnts,1:tail_idx];
        else
            index_out = head_idx:tail_idx;
        end
    end
end

