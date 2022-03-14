function linearity = compute_linearity(pc_scanline,pnt_ref_idx)
    % pc_scanline - (Nx3) points of a single laser scan
    % pnt_ref_ids - index of the point for which smoothness is computed.

    num_pnts  = size(pc_scanline,1);
    pnt_ref   = pc_scanline(pnt_ref_idx,:); % reference point   
%     pc_scanline(pnt_ref_idx,:) = []; % reference is excluded during the computation
%     pnt_ref   = mean(pc_scanline,1);
    
    diff = pc_scanline - pnt_ref;  % calculate the difference between all points in the slide window with the chosen point
    Covar = diff'*diff;
    [~,sv,~] = svd(Covar);

    Covar = Covar./num_pnts;
    % compute singular values & vectors
    [~,sv,~] = svd(Covar);
    eig1 = sv(1,1); 
    eig2 = sv(2,2); 
    eig3 = sv(3,3);

    % Compute Geometric Features (for points in 3D pc)
    linearity  = (eig1 - eig2)/eig1;
%     planarity  = (eig2 - eig3)/eig1;
%     sphericity = eig3/eig1;
%     omnivari   = (eig1*eig2*eig3)^(1/3); % Omnivariance
%     anisotropy = (eig1-eig3)/eig1;
%     eigensum   = eig1+eig2+eig3;
%     surfvar    = eig3/eigensum; % surface variation in 3D
%     eigentropy = -(eig1*log(eig1) + eig2*log(eig2) + eig3*log(eig3)); % Eigen Entropy
    
end