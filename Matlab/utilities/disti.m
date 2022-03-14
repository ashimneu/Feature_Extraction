function dist = disti(idx1,idx2,pc_xyz)
    % pc - [numx3] array of points
    % compute distance of points in pc array.
    % points are selected using their index in the array.
    if nargin == 2
        % dir = '.\data';
        % fname = 'pc_xyz.mat';
        path = 'pc_xyz.mat';
        pc = load(path);
        pc_xyz = pc.pc_xyz;
    end
    pnt1 = pc_xyz(idx1,:);
    pnt2 = pc_xyz(idx2,:);
    dist = norm(pnt2 - pnt1,2); % euclidean distance
end

