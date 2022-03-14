function yhat = proposedfunc(b,x)
    N = size(x,1);
    yhat = zeros(N,1);
    b1 = b(1);
    b2 = b(2);
    b3 = b(3);
    for j = 1:N
        x1 = x(j,1);
        yhat(j) = b3 - exp(-(b1.*x1 + b2).^2);
%         yhat(j) = exp(b1.*x1 + b2*x2 + b3)^2;
    end 
%     yhat = b1 + b2.*x1 + b3.*x2;
end

% function yhat = myfunc3(m,x)
%     b1 = m(1);
%     b2 = m(2);
%     b3 = m(3); 
%     yhat = b1 + b2.*x(:,1) + b3.*x(:,2);
% end


% function yhat = myfunc3(m,x)
%     N = size(x,1);
%     yhat = zeros(N,1);
%     b1 = m(1);
%     b2 = m(2);
%     b3 = m(3);
% 
%     for j = 1:N
%         x1 = x(j,1);
%         x2 = x(j,2);
%         yhat(j) = b1 + b2.*x1 + b3.*x2;
%     end 
% %     yhat = b1 + b2.*x1 + b3.*x2;
% end