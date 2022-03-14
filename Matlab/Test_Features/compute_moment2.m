function [FOFA,FOSA,SOFA,SOSA] = compute_moment2(pc,pnt_ref)
    % For a select point p and its corresponding neighborhood, 
    % various moments are computed.
    % INPUT
    % ref_pnt   - (1x3) reference point
    % pc        - (Nx3) point cloud of neighbors of refernce point
    % OUTPUT
    % FOFA - 1st order, 1st axis,   FOSA - 1st order, 2nd axis
    % SOFA - 2nd order, 1st axis,   SOSA - 2nd order, 2nd axis
    %
    % Following computations are based on moment definitions given in
    % Section 3.2 Feature Extraction & Table 1 in [1].
    % [1] - Fast semantic segmentation of 3D point clouds with strongly
    % varying density by Timo Hackel, Jan Wegner, Konrad Schindler
        
    
    num_pnts  = size(pc,1);
    [med, m_idx] = mfames(pc); % medoid
    diff  = pc - med;
    Covar = diff'*diff./num_pnts;
    
%     % compute eigen values & eigen vectors
%     [Evec, Eval, ~] = eig(Covar);
%     [Eval,Idx]= sort(diag(Eval),'descend');
%     eig1 = Eval(1); evec1 = Evec(:,Idx(1));
%     eig2 = Eval(2); evec2 = Evec(:,Idx(2));
%     eig3 = Eval(3); evec3 = Evec(:,Idx(3));

    % compute singular values & vectors
    [U,S,~] = svd(Covar);
    eig1 = S(1,1); evec1 = U(:,1);
    eig2 = S(2,2); evec2 = U(:,2);
    eig3 = S(3,3); evec3 = U(:,3);

    diff = pc - pnt_ref;
    FOFA_ith = zeros(num_pnts,1); % FOFA ith term
    FOSA_ith = zeros(num_pnts,1); % FOSA ith term
    SOFA_ith = zeros(num_pnts,1); % SOFA ith term
    SOSA_ith = zeros(num_pnts,1); % SOSA ith term
    for idx = 1:num_pnts
        dirvec_i = diff(idx,:)'; % direction vector from reference point along ith point. 
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

