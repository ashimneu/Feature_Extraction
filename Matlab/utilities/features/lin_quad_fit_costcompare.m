function [flag,cost_1storderfit,cost_2ndorderfit] = lin_quad_fit_costcompare(x,y,threshold)
    % flag - 1 - linear fit, 2 - nonlinear fit
    % check linear & quadratic fit over eigensum plot
    if nargin <3
        threshold = 0.1; % default for linear fit detection in SHU pointcloud
    end

    coeff_1st_order = polyfit(x,y,1);
    f1 = polyval(coeff_1st_order,x);
    cost_1storderfit  = sum(abs(y - f1)/numel(f1));
    coeff_2nd_order = polyfit(x,y,2);
    f2 = polyval(coeff_2nd_order,x);
    cost_2ndorderfit  = sum(abs(y - f2)/numel(f2));

    cost_diff = cost_1storderfit - cost_2ndorderfit;

    if cost_diff < threshold
        flag = 1;
    elseif cost_diff >= threshold
        flag = 2;
    end
end

