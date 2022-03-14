function [sorted_vals,sorting_idx] = sort_index(vals)
    % vals - index values. It's entries are indices of
    % points in a point cloud scanline. This function is required when it is
    % found that a chosen cluster consists of points that are from the 
    % beginning & end of a scanline (based on their index in the scanline). 
    % A line plot of index value of such points has a drop from 
    % index=max_index to index=1 which can have issue when processing
    % features and it depends on indices to be serially increasing.
    % To make index values increase sequentially across the pointcloud, 
    % index values in ivals array are adjusted so that index values continue 
    % to increase after index=max_index. For instance, index=1,2,3,.. are 
    % assigned index = max_index+1,max_index+2,max_index+3,...
    
    % val_S = lowest value in shifted_ivals (not first entry)
    % val_L = largest value in shifted_ivals (not last entry)
    
    % sort the values in ascending order
    [sorted_vals,sorting_idx] = sort(vals,'ascend');
    
    % check for large jump in sequentially sorted index values.
    diff_val = sorted_vals(2:end) - sorted_vals(1:end-1);
    [maxdiff,i_maxdiff] = max(diff_val);   
    if maxdiff >= 100
        num = numel(vals); % number of values
        shiftcount = num-i_maxdiff; % count of entries to shift
        sorted_vals = circshift(sorted_vals,shiftcount);
        sorting_idx = circshift(sorting_idx,shiftcount);
    end
end

