function [coeff,res2_sum] = linearfit(xypoints,flags)
    
    pc_xy = xypoints(flags,:);
    N = size(pc_xy,1);
    y = pc_xy(:,2);
    x = pc_xy(:,1);
    H = [x, ones(N,1)];
    coeff = inv(H'*H)*H'*y; % estimated coefficients
    
    yhat = H*coeff;

    residual = y - yhat;
    res2 = residual.^2; % squared residuals
    res2_sum = sum(res2); % sum of squared residuals
end

