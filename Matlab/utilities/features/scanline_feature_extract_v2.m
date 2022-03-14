function [cor_pc,ln_pc,cor_pc_idx,ln_pc_idx,otherfeatures] = scanline_feature_extract_v2(p,pc)
    % p - struct of parameters
    % PointCloud - [numx3] array of points obtained from one laser scanline
    % ndef - neighbor definition. accepted args - 'knn','inradius'
    % cor_pc - corner points
    % ln_pc  - linear points
    % otherfeatures - compute value of various geomteric features
    % Cell in otherfeatures are sequentially indexed to denote
    % L,P,S,O,A,Sum,SV,Entropy,FOFA,FOSA,SOFA,SOSA, respectively.
    % L - Linearity,     P - Planarity,  S - Sphericity
    % O - Ominivariance, A - Anisotropy, Sum - sum of eigen vals
    % SV - Surface Variation, Entropy - Eigen Entropy
    % FOFA - 1st order, 1st axis,   FOSA - 1st order, 2nd axis
    % SOFA - 2nd order, 1st axis,   SOSA - 2nd order, 2nd axis
    
    %     if lower(ndef) == 'knn'
    %     
    %     elseif lower(ndef) == 'inradius'
    %         [indices,dists] = findNeighborsInRadius(ptxyz,point,radius);
    %     else
    %         disp('Value passed in ''ndef'' is invalid.')
    %         disp('Accepted entries for ndef are ''knn'' & ''inradius''.')
    %         error('')
    %     end
    
    window_def = p.sw_def; % 1,2,3 are knnSearch v1,v2 and distance based window respectively

    valid_range = inf; % threshold of valid slide window
    sz_window = p.sz_window; % # of points in a sliding window ( should be odd integer)
    hf_sz = (sz_window-1)/2;
      
    num_pnts = size(pc,1);
    if(sz_window>num_pnts)
        disp('window too larger or equal to the number of points')
    end
    
    % empty variables to store computed feature values
    feature_pc = [pc, zeros(num_pnts,1),zeros(num_pnts,1)]; % default type as 0;
    [L,P,S,O,A,Sum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dp,len_pr,len_sc] = deal(nan(num_pnts,1));
    otherfeatures = cell(1,16); % number of features to extract
    
    sz_pc_sw = zeros(num_pnts,1);
    for i = 1:num_pnts
        switch window_def
            case 1
                pc_sw = SWindowknnSearch_v1(pc,i,hf_sz);
                sz_pc_sw(i) = sz_window;
            case 2
                pc_sw = SWindowKnnSearch_v2(pc,i,hf_sz);
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
                
                % compute covariance matrix
                diff = pc_sw - medoid;
                H = diff'*diff;
                Covar = H./sz_pc_sw(i);

                % compute singular values & vectors of covar. matrix
                [~,sv,~] = svd(Covar);
                eig1 = sv(1,1); 
                eig2 = sv(2,2); 
                eig3 = sv(3,3);
        
                % Compute Geometric Features (for points in 3D pc)
                L(i)   = (eig1 - eig2)/eig1;
                P(i)   = (eig2 - eig3)/eig1;
                S(i)   = eig3/eig1;
                O(i)   = (eig1*eig2*eig3)^(1/3);
                A(i)   = (eig1-eig3)/eig1;
                Sum(i) = eig1+eig2+eig3;
                SV (i) = eig3/Sum(i); % surface variation in 3D
                Entropy(i) = -(eig1*log(eig1) + eig2*log(eig2) + eig3*log(eig3)); % Eigen Entropy
                
                % compute moments
                [FOFA(i), FOSA(i), SOFA(i), SOSA(i)] = compute_moment2(pc_sw,pnt_ref);
                
                % compute smoothness
                Sm(i) = compute_smoothness(pc_sw,pnt_ref_idx);
                
                % compute dot product by fitting st. line on either side of
                % current point
                [dp(i),~,len_pr(i),len_sc(i)] = dotprod_fitlines(p,pc_sw,pnt_ref_idx);

                if p.dcsnfeat_cor == 1
                    % linearity to decide if ith point is a corner or
                    % not corner
                    threshold_L = p.threshold_L; % Linearity threshold
                    if L(i) <= threshold_L
                        feature_pc(i,5) = 2; % corner
                    elseif L(i) > threshold_L
                        feature_pc(i,5) = 1; % not corner
                    end
                elseif p.dcsnfeat_cor == 2
                    % smoothness to decide if ith point is a corner or
                    % not corner
                    threshold_Sm = p.threshold_Sm; % Smoothness threshold
                    if Sm(i) <= threshold_Sm
                        feature_pc(i,5) = 1; % not corner
                    elseif Sm(i) > threshold_Sm
                        feature_pc(i,5) = 2; % corner
                    end
                elseif p.dcsnfeat_cor == 3
                    % use threshold to determine whether a point cloud has a corner
                    % if it only consists of points close to being a line
                    threshold_Sm = p.threshold_Sm; % Smoothness threshold
                    threshold_L  = p.threshold_L;  % Linearity threshold
                    temp1 = 0; 
                    temp2 = 0;
                    if L(i) <= threshold_L
                        temp1 = 2; % corner
                    elseif L(i) > threshold_L
                        temp1 = 1; % not corner
                    end
                    if Sm(i) <= threshold_Sm
                        temp2 = 1; % not corner
                    elseif Sm(i) > threshold_Sm
                        temp2 = 2; % corner
                    end
                     feature_pc(i,5)
                end
            end
        end
        %     [coeff, score, latent] = pca(pc_sw)    
    end
    
    ln_pc_idx  = find(feature_pc(:,5)==1);
    cor_pc_idx = find(feature_pc(:,5)==2);
    ln_pc  = feature_pc(ln_pc_idx,1:3);
    cor_pc = feature_pc(cor_pc_idx,1:3);
    [otherfeatures{:}] = deal(L,P,S,O,A,Sum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dp,len_pr,len_sc);
end