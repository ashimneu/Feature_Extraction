function [features] = compute_feature_using_sw_v2(p,PointCloud,sz_window)
    % PointCloud - [Nx3] set of points from a single laser scan
    % sz_window - size of Sliding Window (SW). Should be an odd integer.
    % features - entries are 'linearity', 'planarity', 
    %               'sphericity', 'omnivar', 'anisotropy','eigensum',
    %               'surfvar', 'eigenentropy', 'smoothness','FOFA',
    %               'FOSA','SOFA','SOSA','dp','len_pr','len_sc'.
    % [1] - FAST SEMANTIC SEGMENTATION OF 3D POINT CLOUDS WITH STRONGLY
    %       VARYING DENSITY by Hackel, Wegner & Schindler
    % [2] - LOAM: Lidar Odometry and Mapping in Real-time by Zhang & Singh
    
    window_def = p.sw_def; % 1,2,3 are knnSearch v1,v2 and distance based window respectively
    EXCLUDE_SW_OF_ENDPOINTS = true;
    
    if mod(sz_window,2) == 0
        error('Value of sz_window needs to be an odd integer.')
    end

    pc = PointCloud;
    num_pnts = size(pc,1);

    % empty variables to store computed feature values
    [L,P,S,O,A,Esum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dp,len_pr,len_sc,angle] = deal(nan(num_pnts,1));
    features = cell(1,16); % number of features to extract
    sz_pc_sw = zeros(num_pnts,1); % SW size for neighborhood of each point

    % default min & max values of SW; odd integers only accepted.
    sz_min_max = [5,odd_int(num_pnts)]; % min & max        
    
    for i = 1:num_pnts
        switch window_def
            case 1
                sz_window = checkWindowSize(sz_window,num_pnts,sz_min_max);   
                hf_sz = (sz_window-1)/2;
                pc_sw = SWindowknnSearch_v1(pc,i,hf_sz);
                sz_pc_sw(i) = sz_window;
            case 2
                sz_window = checkWindowSize(sz_window,num_pnts,sz_min_max);   
                hf_sz = (sz_window-1)/2;
                pc_sw = SWindowknnSearch_v2(pc,i,hf_sz);
                sz_pc_sw(i) = sz_window;
            case 3
                [pc_sw,pnt_ref_i] = SWindowInRadiusSearch(p.radius,pc,i);
                sz_pc_sw(i) = size(pc_sw,1);
        end

        if ~isempty(pc_sw)
            dis_sep = norm(pc_sw(1,:)-pc_sw(end,:));
            valid_range = inf;
            if (dis_sep < valid_range) % only analyze when all points are within a valid small range of slide window
                switch window_def
                case 1
                    pnt_ref = pc(i,:); % reference point USE ME if using SWindowknnSearch_v1/v2 
                    pnt_ref_idx = hf_sz+1; % index of the central point in the sliding window; USE ME if using SWindowknnSearch_v1/v2
                    %[medoid,~] = mfames(pc_sw'); % true medoid                 
                    % using mid point in sliding window as close approximation 
                    % of medoid since medoid calculation is computationally expensive.
                    medoid = pc_sw(hf_sz,:);
                case 2
                    pnt_ref = pc(i,:); % reference point USE ME if using SWindowknnSearch_v1/v2
                    pnt_ref_idx = hf_sz+1; % index of the central point in the sliding window; USE ME if using SWindowknnSearch_v1/v2
                    %[medoid,~] = mfames(pc_sw'); % true medoid 
                    medoid = pc_sw(hf_sz,:); % mid point of sliding window. See case 1 for more details.
                case 3
                    pnt_ref_idx = pnt_ref_i;
                    pnt_ref = pc(pnt_ref_idx,:); % reference point USE ME if using SWindowInRadiusSearch
                    [medoid,~] = mfames(pc_sw'); % medoid (computationally expensive)
                    medoid = medoid'; % convert to row vector
                    % medoid = pc_sw(hf_sz,:);
                end 
                
                % Compute Geometric features defined in [1], table 1.
                [L(i),P(i),S(i),O(i),A(i),Esum(i),SV(i),Entropy(i),Sm(i),FOFA(i),FOSA(i),SOFA(i),SOSA(i)] = compute_features(pc_sw,pnt_ref_idx);
                p.i = i;
                [dp(i),angle(i),len_pr(i),len_sc(i)] = dotprod_fitlines_v1(p,pc_sw,pnt_ref_idx);
                % pnt_ref = pc(i,:);
                % [FOFA(i),FOSA(i),SOFA(i),SOSA(i)] = compute_moment2(pc_sw,pnt_ref);
            end
        end
    end
    [features{:}] = deal(L,P,S,O,A,Esum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dp,len_pr,len_sc);
end

function [L,P,S,O,A,Esum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA] = compute_features(pc_sw,pnt_ref_idx)
    % pc_scanline - (Nx3) points of a single laser scan
    % pnt_ref_idx - index of reference point in pc_sw, for which 
    %               geometric features are computed.
    % [1] - FAST SEMANTIC SEGMENTATION OF 3D POINT CLOUDS WITH STRONGLY
    %       VARYING DENSITY by Hackel, Wegner & Schindler
    % [2] - LOAM: Lidar Odometry and Mapping in Real-time by Zhang & Singh
    
    pnt_ref  = pc_sw(pnt_ref_idx,:); % reference point

    num_pnts = size(pc_sw,1);
    [medoid,~] = mfames(pc_sw'); % medoid
    
    % compute covariance matrix
    diff = pc_sw - medoid';  % difference between pc points with medoid
    Covar = diff'*diff;
    Covar = Covar./num_pnts;

    % compute singular values & vectors
    [U,sv,~] = svd(Covar);
    eig1 = sv(1,1); evec1 = U(:,1);
    eig2 = sv(2,2); evec2 = U(:,2);
    eig3 = sv(3,3);
    
    % Compute Geometric features defined in [1], table 1.
    L  = (eig1 - eig2)/eig1; % Linearity
    P  = (eig2 - eig3)/eig1; % Planarity
    S  = eig3/eig1;          % Sphericity
    O  = (eig1*eig2*eig3)^(1/3); % Omnivariance
    A  = (eig1-eig3)/eig1;   % Anisotropy
    Sm = compute_smoothness(pc_sw,pnt_ref_idx); % Smoothness (defined in [2],eqn(1))
    Esum = eig1+eig2+eig3;   % Sum of Eigen Values
    SV = eig3/Esum;          % surface variation
    Entropy = -(eig1*log(eig1) + eig2*log(eig2) + eig3*log(eig3)); % Eigen Entropy    
    
    % compute moment features    
    diff2    = pc_sw - pnt_ref; % difference between pc points with reference point
    FOFA_ith = zeros(num_pnts,1); % FOFA ith term
    FOSA_ith = zeros(num_pnts,1); % FOSA ith term
    SOFA_ith = zeros(num_pnts,1); % SOFA ith term
    SOSA_ith = zeros(num_pnts,1); % SOSA ith term    
    for idx = 1:num_pnts
        dirvec_i = diff2(idx,:)'; % direction vector from reference point along ith point. 
        FOFA_ith(idx) = dot(dirvec_i,evec1); 
        FOSA_ith(idx) = dot(dirvec_i,evec2);  
        SOFA_ith(idx) = dot(dirvec_i,evec1)^2;
        SOSA_ith(idx) = dot(dirvec_i,evec2)^2;
    end
    FOFA = sum(FOFA_ith); % 1st order, 1st axis
    FOSA = sum(FOSA_ith); % 1st order, 2nd axis
    SOFA = sum(SOFA_ith); % 2nd order, 1st axis
    SOSA = sum(SOSA_ith); % 2nd order, 2nd axis
end
