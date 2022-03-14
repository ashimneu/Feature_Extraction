function [FOFA,FOSA,SOFA,SOSA]  = moments(pc,idx_pnt_ref,U)
    % For a select point p and its corresponding neighborhood, 
    % various moments are computed.
    % INPUT:
    % pnt_ref - index of reference point in the pc array
    % pc      - (Nx3) neighboring points of reference point
    % U       - Matrix of orthogonal eigen vectors of Covariance matrix
    % OUTPUT:
    % FOFA - 1st order, 1st axis,   FOSA - 1st order, 2nd axis
    % SOFA - 2nd order, 1st axis,   SOSA - 2nd order, 2nd axis
    %
    % Following computations are based on moment definitions given in
    % Section 3.2 Feature Extraction & Table 1 in [1].
    % [1] - Fast semantic segmentation of 3D point clouds with strongly
    % varying density by Timo Hackel, Jan Wegner, Konrad Schindler
        
    
    num_pnts  = size(pc,1);
    pnt_ref = pc(idx_pnt_ref,:);
    evec1 = U(:,1);
    evec2 = U(:,2);
    evec3 = U(:,3);

    % compute moment features    
    diff2    = pc - pnt_ref; % difference between pc points with reference point
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

