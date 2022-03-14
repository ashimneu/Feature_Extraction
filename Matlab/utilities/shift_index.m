function [shifted_ivals] = shift_index(ivals)
    % ivals is an acronym for index values. It's entries are indices of
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
    
    diff_ivals = ivals(2:end) - ivals(1:end-1);
    [~,i_maxdiff] = max(diff_ivals);
    max_ivals = max(ivals);
    val_L = max_ivals + i_maxdiff;
    new_vals = max_ivals+1:1:val_L;
    shifted_ivals = ivals;
    shifted_ivals(1:i_maxdiff) = new_vals;

%     [val_S,loc_S] = min(ivals); % location & value of smallest entry in ivals array
%     [val_L,loc_L] = max(ivals); % location & value of largest entry in ivals array
%     num  = numel(ivals); % total number of entries
%     diff = num - loc_L; % number of entries to be replaced with updated value in ivals array.
%     shifted_ivals = ivals;
%     shifted_ivals(loc_S:end) = val_L + 1:1:diff;
end

