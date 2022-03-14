function O = omnivariance(sv)
    % sv - Matrix of singular values of a covariance matrix
    eig1 = sv(1,1);
    eig2 = sv(2,2);
    eig3 = sv(3,3);
    
    O  = (eig1*eig2*eig3)^(1/3); % Omnivariance
end

