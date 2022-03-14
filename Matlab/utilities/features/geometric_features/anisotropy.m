function A = anisotropy(sv)
    % sv - Matrix of singular values of a covariance matrix
    eig1 = sv(1,1);
    eig2 = sv(2,2);
    eig3 = sv(3,3);
    
    A  = (eig1-eig3)/eig1; % Anisotropy
end

