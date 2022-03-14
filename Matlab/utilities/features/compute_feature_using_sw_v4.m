function [f] = compute_feature_using_sw_v4(p,PointCloud,sz_window)
    % PointCloud - [Nx3] set of points from a single laser scan
    % sz_window - size of Sliding Window (SW). Should be an odd integer.
    % f - struct whose members are 'linearity', 'planarity', 
    %               'sphericity', 'omnivar', 'anisotropy','eigensum',
    %               'surface variation', 'eigenentropy', 'smoothness',
    %               'FOFA', 'FOSA','SOFA','SOSA',
    %               'dp','angle','len_pr','len_sc'.
    % [1] - FAST SEMANTIC SEGMENTATION OF 3D POINT CLOUDS WITH STRONGLY
    %       VARYING DENSITY by Hackel, Wegner & Schindler
    % [2] - LOAM: Lidar Odometry and Mapping in Real-time by Zhang & Singh
    
    window_def = p.sw_def; % 1,2,3 are knnSearch v1,v2 and distance based window respectively
    
    if mod(sz_window,2) == 0
        error('Value of sz_window needs to be an odd integer.')
    end

    pc = PointCloud;
    num_pnts = size(pc,1);

    % declare variables to store computed feature values
    [f.L,f.P,f.S,f.O,f.A,f.Esum,...
    f.SV,f.Entropy,f.Sm,...
    f.FOFA,f.FOSA,f.SOFA,f.SOSA,...
    f.dotprod,f.angle,f.len_pr,f.len_sc] = deal(nan(num_pnts,1));

    sz_pc_sw = zeros(num_pnts,1); % SW size for neighborhood of each point

    % default min & max values of SW; odd integers only accepted.
    sz_min_max = [5,odd_int(num_pnts)]; % min & max        
    
    for i = 1:num_pnts
        p.idx_debug = p.idxa(i); % DELETE ME AFTER DEBUGGING!
        
        switch window_def
            case 1
                sz_window = checkWindowSize(sz_window,num_pnts,sz_min_max);   
                hf_sz = (sz_window-1)/2;
                [pc_sw,idx_pnt_ref_sw,idx_out] = SWindowknnSearch_v1(pc,i,hf_sz);
                sz_pc_sw(i) = sz_window; % number of points in SW
                p.idxa_sw = p.idxa(idx_out);
            case 2
                sz_window = checkWindowSize(sz_window,num_pnts,sz_min_max);   
                hf_sz = (sz_window-1)/2;
                [pc_sw,idx_pnt_ref_sw,idx_out] = SWindowknnSearch_v2(pc,i,hf_sz);
                sz_pc_sw(i) = sz_window; % number of points in SW
                p.idxa_sw = p.idxa(idx_out);
            case 3
                [pc_sw,idx_pnt_ref_sw,idx_out] = SWindowInRadiusSearch(p.radius,pc,i); % idx_pnt_ref_sw - index of refernce in pc_sw array
                sz_pc_sw(i) = size(pc_sw,1); % number of points in SW
                p.idxa_sw = p.idxa(idx_out);
        end

        if ~isempty(pc_sw)
            dis_sep = norm(pc_sw(1,:)-pc_sw(end,:));
            valid_range = inf;
            if (dis_sep < valid_range) % only analyze when all points are within a valid small range of slide window
                switch window_def
                case 1
                    pnt_ref_idx = idx_pnt_ref_sw; % index of the central point in the sliding window; USE ME if using SWindowknnSearch_v1/v2
                    %[medoid,~] = mfames(pc_sw'); % true medoid                 
                    % using mid point in sliding window as close approximation 
                    % of medoid since medoid calculation is computationally expensive.
                    medoid = pc_sw(hf_sz,:);
                case 2
                    pnt_ref_idx = idx_pnt_ref_sw; % index of refence point in the SW
                    %[medoid,~] = mfames(pc_sw'); % true medoid 
                    medoid = pc_sw(hf_sz,:); % mid point of sliding window. See case 1 for more details.
                case 3
                    
                    pnt_ref_idx = idx_pnt_ref_sw;                    
                    % uncomment this if true medoid is preferred for feature computation.
                    [medoid,~] = mfames(pc_sw'); % medoid (computationally expensive)
                    medoid = medoid'; % convert to row vector
                    
                    % uncomment this if approximation of mediod is preferred for feature computation.
                    % medoid = pc_sw(floor(sz_pc_sw(i)/2),:); % rough approximation of medoid 
                end 
               
                % compute covariance of points in SW
                diff = pc_sw - medoid;  % difference between pc points with medoid
                Covar = diff'*diff;
                Covar = Covar./sz_pc_sw(i);
            
                % compute singular values & vectors
                [U,sv,~] = svd(Covar);

                % Compute Geometric features defined in [1], table 1.
                f.L(i) = linearity(sv);
                f.P(i) = planarity(sv);
                f.S(i) = sphericity(sv);
                f.O(i) = omnivariance(sv);
                f.A(i) = anisotropy(sv);
                f.Esum(i) = eigen_sum(sv);
                f.SV(i) = surfacevariance(sv);
                f.Entropy(i) = eigen_entropy(sv);
                f.Sm(i) = smoothness(pc_sw,pnt_ref_idx);
                [f.FOFA(i),f.FOSA(i),f.SOFA(i),f.SOSA(i)] = moments(pc_sw,pnt_ref_idx,U);

                % compute dot product, length of fitted lines & angle between them.
                [f.dotprod(i),f.angle(i),f.len_pr(i),f.len_sc(i)] = dotprod_fitlines_v3(p,pc_sw,pnt_ref_idx);                
            end
        end
    end
end