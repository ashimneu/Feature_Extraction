function [FOFA,FOSA,SOFA,SOSA] = compute_moment_using_sw(PointCloud,sz_window)
    % Following moments of various order along various axes (denoted by eig vector 
    % correponding to 1st & 2nd largest eig. vals.) are computed for a reference point p 
    % in its neighborhood (determined by a sliding window).
    
    % INPUT
    % PointCloud - (Nx3) point cloud of 1 laser taken as single scanline.
    %
    % OUTPUT
    % FOFA - 1st order, 1st axis
    % FOSA - 1st order, 2nd axis
    % SOFA - 2nd order, 1st axis
    % SOSA - 2nd order, 2nd axis
    %
    % Sliding window (sw) of odd integer length traverses along the
    % scanline and select points to define neighborhood of a reference point.
    % sliding window begins   
    
    % since sliding window is centered on a reference point and spans
    % equally (up to window_size/2 number of points) on both sides of the 
    % PointCloud array. For window_size/2 number of points on leading &
    % trailing end of a scanline, getting such sliding window is not
    % feasible. Hence, sliding window based neighboorhood point extraction
    % is avoided for these points on leading & trailing end of the scanline
    % given in PointCloud.
   
    % at begining & end of scanline, sliding window cannot be centered
    % about the points 

    EXCLUDE_SW_OF_ENDPOINTS = true;

    % pc =  PointCloud.Location;
    pc = PointCloud;
    num_pnts = size(pc,1);

    % declare variables to store the computed values.
    FOFA = nan(num_pnts,1);
    FOSA = nan(num_pnts,1);
    SOFA = nan(num_pnts,1);
    SOSA = nan(num_pnts,1);
    
%     sz_window = 19; % Size of the Sliding Window (SW) of each make sure to be odd integer
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
                pnt_ref = pc(i,:);
                [FOFA(i), FOSA(i), SOFA(i), SOSA(i)] = compute_moment2(pc_sw,pnt_ref);
            end
        end
    end    
end