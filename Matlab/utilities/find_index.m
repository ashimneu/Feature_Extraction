function [idx_s,idx_g,pc_s,pc_g] = find_index(feature_value,threshold,pc)
    % get indices of points whose 'feature_value's are either greater or less
    % compared to given threshold.
    idx_s = find(feature_value <= threshold); 
    idx_g = find(feature_value > threshold);
    pc_s = pc(idx_s,:); % points whose feature value <= threshold
    pc_g = pc(idx_g,:); % points whose feature value <= threshold
end