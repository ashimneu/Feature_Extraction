function [smoothness] = compute_smoothness_using_sw(PointCloud,sz_window)
    % PointCloud - [Nx3] set of points from a single laser scan
    % sz_window - size of Sliding Window (SW). Should be an odd integer.
    if mod(sz_window,2) == 0
        error('Value of sz_window needs to be an odd integer.')
    end

    EXCLUDE_SW_OF_ENDPOINTS = true;

    % pc =  PointCloud.Location;
    pc = PointCloud;
    num_pnts = size(pc,1);

    % declare variables to store the computed values.
    smoothness = nan(num_pnts,1);
    
    % default min & max values of SW.
    sz_window_min = 5;
    sz_window_max = num_pnts;
    hf_sz = (sz_window-1)/2;    

    if (sz_window > num_pnts)
        disp('Sliding window size is larger than the number of points.')
        disp('Setting window size to max. number of points.')
        sz_window = sz_window_max;
    elseif (sz_window < sz_window_min)
        disp('Sliding window size is less than threshold of 5.')
        disp('Setting window size to 5.')
        sz_window = sz_window_min;
    end
    
    for i = 1:num_pnts
        sw_start_idx = i - hf_sz; % index of initial point in sliding window
        sw_stop_idx  = i + hf_sz; % index of last point in sliding window
        pc_sw = [];
        if (sw_start_idx <= 0) && ~EXCLUDE_SW_OF_ENDPOINTS % check if start index is at leading end of pointcloud scanline
            sw_start_idx = 1;
            sw_stop_idx  = sz_window;
            pc_sw = pc(sw_start_idx:sw_stop_idx,:); 
        elseif (sw_stop_idx > num_pnts) && ~EXCLUDE_SW_OF_ENDPOINTS % check if stop index is at training portion of pointcloud scanline
            sw_start_idx = num_pnts - sz_window + 1;
            sw_stop_idx  = num_pnts;
            pc_sw = pc(sw_start_idx:sw_stop_idx,:); 
        elseif  (sw_start_idx > 0) && (sw_stop_idx <= num_pnts)
            sw_start_idx = i - hf_sz;
            sw_stop_idx  = i + hf_sz;
            pc_sw = pc(sw_start_idx:sw_stop_idx,:); 
        end           
        
        if ~isempty(pc_sw)
            dis_sep = norm(pc_sw(1,:)-pc_sw(end,:));
            valid_range = inf;
            if (dis_sep < valid_range) % only analyze when all points are within a valid small range of slide window
                pnt_ref_idx = hf_sz+1; % index of the central point in the sliding window
                smoothness(i) = compute_smoothness(pc_sw,pnt_ref_idx);
            end
        end
    end    
end
