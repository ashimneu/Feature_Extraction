function L = linearity(sv)
    % sv - Matrix of singular values of a covariance matrix
    eig1 = sv(1,1);
    eig2 = sv(2,2);
    
    L  = (eig1 - eig2)/eig1; % Linearity
end

