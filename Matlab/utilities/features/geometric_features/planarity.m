function P = planarity(sv)
    % sv - Matrix of singular values of a covariance matrix
    eig1 = sv(1,1);
    eig2 = sv(2,2);
    eig3 = sv(3,3);
    
    P  = (eig2 - eig3)/eig1; % Planarity
end

