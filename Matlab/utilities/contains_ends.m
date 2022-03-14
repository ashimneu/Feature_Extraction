function output = contains_ends(indices,index_endpoint)
    % output - 1 or 0
    
    flag_idx_initpoint = indices == 1; % checking for index of initial point
    contains_initpoint = sum(flag_idx_initpoint);   % 1 if index is found otherwise 0
    flag_idx_endpoint  = indices == index_endpoint; % checking for index of final
    contains_endpoint  = sum(flag_idx_endpoint);    % 1 if index is found otherwise 0
    
    if contains_initpoint && contains_endpoint
        % if both inital & final points' index are present.
        output = 1;
    else
        % if either of inital or final points' indices aren't present.
        output = 0;
    end    
end

