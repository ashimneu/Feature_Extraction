function [density,mean_dist,dist] = compute_density(pc,n)
    % pc - [Nx3] set of points from a single laser scan
    % n  - number of points over which local density is computed.
    %      n(>=3) should be odd (n = 2m + 1).
    % density - at a chosen point, it is computed by picking equal number
    %           of preceeding & succeeding points
    % dist - distance between point i & point i+1;
    % mean_distance - average of all dist values (possibly trimmed)
       
    num_pnts = size(pc,1);
    
    % compute euclidean distance between succesive points
    distance = zeros(num_pnts-1,1);
    for i = 1:num_pnts-1
        distance(i) = norm(pc(i+1,:) - pc(i,:),2);
    end
    
    % compute local density, given minimum number of points    
    m = (n-1)/2;
    % compute density = # points per meter
    density = nan(num_pnts,1);
    i_start = 1 + m;
    i_end   = num_pnts - m;
    if (i_end - i_start) >= 1
        for i = i_start:i_end
            start_pnt_idx = i - m;
            end_pnt_idx   = i + m;
            local_line_length = sum(distance(start_pnt_idx:end_pnt_idx-1));
            density(i) = n/local_line_length;
        end
    end

    if nargout > 1
        % mean distance between points in pc
        trim_percent = 5;
        mean_dist = trimmean(distance,trim_percent); % trimmed for robustness against some outliers
        dist = distance; 
        dist(end+1) = nan;
    end
end