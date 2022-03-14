function [medoid medoid_index] = mfames(cluster_members, ...
                                        num_axes, ...
                                        num_candidates, ...
                                        dist_fn)
%%%
% An implementation of the (Mutliple) FAst MEdoid Selection
% algorithm (MFAMES) [Paterlini et al. 2011]
%
% This algorithm is nondeterministic, so be sure to seed the PRNG
% appropriately.
%
% Author: Stephen Pratt (sdp2128@gmail.com)
%
% Args:
%   cluster_members: Ideally a 1 x N cell array where
%   cluster_member{1, i} is a column vector representing the i'th
%   member of the cluster whose medoid is sought. Also accepts a
%   d x N matrix (which will be converted into a cell array).
%
%   num_axes (optional): Number of pivot axes to consider during 
%   computation. Useful for high-dimensional, irregular, or
%   spiculated datasets. Efficiency scales linearly with this
%   value. Defaults to 2.
%
%   num_candidates (optional): Number of medoid candidates to
%   return. If num_candidates is 3, then an estimate of the 3 best
%   medoid candidates will be returned in decending order of
%   confidence. Defaults to 1.
%
%   dist_fn (optional): Function to use as a dissimilarity
%   metric. Should be a handle to a function that accepts two
%   column vectors as input and returns a real scalar distance
%   measure as output. Defaults to @(a, b)(norm(a - b))
%
%
% Returns:
%   The index of the medoid column in cluster_members
%%%

if nargin < 2
    num_axes = 2;
end
if nargin < 3
    num_candidates = 1;
end
if nargin < 4
    dist_fn = @(a, b)(norm(a -b));
end

[d N] = size(cluster_members);

if num_candidates > N
    error('num_candidates cannot be larger than size of dataset');
end

if ~iscell(cluster_members)
    cluster_members = mat2cell(cluster_members, d, ones(1, N));
else
    d = size(cluster_members{1}, 1);
end

s_z = cluster_members{randsample(N, 1)};

% s_prime(:,i) and s_dprime(:,i) are the endpoints of the i'th axis
s_prime = NaN(d, num_axes);
s_dprime = NaN(d, num_axes);

% m(i) will be the median distance from s_prime(:,i) that can be
% created by projecting each of the cluster members onto the axis
% created by s_prime(:,i) and s_dprime(:,i)
m = NaN(1, num_axes);

% chord will become the distance between s_prime(1) and s_dprime(1)
chord = NaN;

% Build our axes
for axis_num = 1:num_axes

    if axis_num == 1
        % Compute the endpoints of the first axis

        % Find the point furthest from our randomly selected
        % starter point. This is our first pivot
        [max_dist argmax_dist] = max(cellfun(...
            @(s)(dist_fn(s_z, s)), ...
            cluster_members));
        s_prime(:,1) = cluster_members{argmax_dist};

        % Our second pivot is the point furthest from the first.
        [max_dist argmax_dist] = max(cellfun(...
            @(s)(dist_fn(s_prime(:,1), s)), ...
            cluster_members));
        s_dprime(:,1) = cluster_members{argmax_dist};
        chord = dist_fn(s_prime(:,1), s_dprime(:,1));
    else
        % Use earlier axes to compute the succeeding endpoints

        % pivot_proximity gives an estimate of how far a datapoint
        % is from a pair of pivots
        pivot_proximity = @(s, s_p, s_dp)(...
            abs(chord - dist_fn(s_p, s)) + ...
            abs(chord - dist_fn(s_dp, s)));

        % sum_pivot_proximity sums pivot_proximity for all known pivots
        sum_pivot_proximity = @(s)(sum(...
            arrayfun(@(i)(...
                pivot_proximity(s, s_prime(:,i), s_dprime(:,i))), ...
                1:(axis_num-1))));

        % Our new pivot should be the maximum possible distance
        % from all other pivots, so we minimize sum_pivot proximity.
        [max_dist argmax_dist] = min(cellfun(sum_pivot_proximity, ...
                                             cluster_members));
        s_prime(:,axis_num) = cluster_members{argmax_dist};

        % Now we need to include the newly selected pivot in our
        % sum_pivot_proximity meaure.
        sum_pivot_proximity = @(s)(...
            sum_pivot_proximity(s) + ...
            abs(chord - dist_fn(s_prime(:,axis_num), s)));

        [max_dist argmax_dist] = min(cellfun(sum_pivot_proximity, ...
                                             cluster_members));
        s_dprime(:,axis_num) = cluster_members{argmax_dist};
    end

    % Finally, we compute our projection point on the newly
    % discovered axis.
    s_p = s_prime(:,axis_num);
    s_dp = s_dprime(:,axis_num);

    % Project all datapoints onto the axis and find the median
    % projected distance from the axis endpoints
    proj_dist_fn = @(s)(...
        (dist_fn(s_p, s)^2 + dist_fn(s_p, s_dp)^2 - ...
         dist_fn(s_dp, s)^2) / (2*dist_fn(s_p, s_dp)));

    proj_dists = sort(cellfun(proj_dist_fn, cluster_members));
    m(axis_num) = proj_dists(ceil(N/2));
end


% Sum of differences function describes distance of cluster points
% to a center pivot point (created implicitly by finding our median
% projected distance). s is the cluster point to test, s_p and s_dp
% are axis endpoints, and m is the median projected distance found.
sod_fn = @(s, s_p, s_dp, m)(...
    abs(dist_fn(s_p, s) - m) + ...
    abs(dist_fn(s_dp, s) - (dist_fn(s_p, s_dp) - m)));

% Sigma sum of diffences function is a gives a summation of the sum
% of differences function over all axis values.
sigma_sod_fn = @(s)(sum(...
    arrayfun(@(i)(sod_fn(s, s_prime(:,i), s_dprime(:,i), m(i))), ...
             1:num_axes)));

if num_candidates == 1
    % medoid estimate is the point that minimizes the sigma sum of
    % differences function 
    [min_dist medoid_index] = min(cellfun(sigma_sod_fn, ...
                                          cluster_members));
    medoid = cluster_members{medoid_index};
else
    % select minimum num_candidates values
    [min_dists candidates] = sort(cellfun(sigma_sod_fn,...
                                          cluster_members));
    medoid_index = candidates(:, 1:num_candidates);
    medoid = cluster_members(medoid_index);
end
