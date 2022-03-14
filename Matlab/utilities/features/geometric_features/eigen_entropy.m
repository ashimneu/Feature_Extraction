function entropy = eigen_entropy(sv)
    % sv - Matrix of singular values of a covariance matrix
    eig1 = sv(1,1); 
    eig2 = sv(2,2);
    eig3 = sv(3,3);
    
    entropy = -(eig1*log(eig1) + eig2*log(eig2) + eig3*log(eig3)); % Eigen Entropy    
    
end

