function yhat = proposedfuncnoisy(b,x,std)
    N = size(x,1);
    yhat = zeros(N,1);
    b1 = b(1);
    b2 = b(2);
    b3 = b(3);
    for j = 1:N
        x1 = x(j,1);
        noise = randn*std;
        yhat(j) = b3 - exp(-(b1.*x1 + b2).^2) + noise;
    end 
end
