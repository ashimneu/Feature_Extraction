function [length] = computelength(pc)
    % computes distance between successive points and adds them to yield
    % the sum as the length of given pc.
    % points need to be in series

    num_pnt = size(pc,1);
    diffvec = pc(2:end,:) - pc(1:end-1,:);
    dist = zeros(num_pnt-1,1);
    for j = 1:num_pnt-1
        dist(j) = norm(diffvec(j),2);
    end
    length = sum(dist);
end

