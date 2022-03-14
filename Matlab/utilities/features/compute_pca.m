function [V,E,VT] = compute_pca(pc)
    % pc - [3xnum] set of points
    % V  - Matrix whose columns are orthogonal singular vectors of
    %     Covariance of pc
    % E  - Matrix of singular values of Covariance of pc.
    % VT - Transpose of V 
    
    num  = size(pc,2);
    cent = mean(pc,2); % centroid
    diff = pc - cent;
    Cov  = (diff*diff')./(num-1); % Covariance of pc
    [V,E,VT] = svd(Cov);
end

