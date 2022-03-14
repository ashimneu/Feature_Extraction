function [med,idx] = medoid(points)
    % points -[num x 3] set of points    
    % Author- Ashim Neupane, email: aneup001@ucr.edu
    % Remark: For quicker medoid computation, use MFAMES by Stephen Pratt
    %         instead of this function. MFAMES is available in the file
    %         exchange library of MATLAB Central.

    num = size(points,1); % number of points in the pointcloud   
    norm_2 = @(p,q) norm(p-q); % compute 2-norm (i.e. distance)
    Dist = zeros(num,num);
    for row = 1:num
        for col = row:num
            Dist(row,col) = norm_2(points(row,:),points(col,:));
        end
    end
    Dist = triu(Dist)' + Dist; % make the matrix symmetric
    TotalDist = sum(Dist,1);
    [~,idx_min_dist] = min(TotalDist);
    
    med = points(idx_min_dist,:); % medoid
    idx = idx_min_dist ; % index of medoid in the array of points
end
