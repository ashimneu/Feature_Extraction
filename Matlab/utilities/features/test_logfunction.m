% clear; clc;

% generate raw data
inc_x = 0.001;
data_x = inc_x:inc_x:100;
true_a = 0.05;
true_b = -0.5;
data_y = -exp(-(true_a.*data_x + true_b).^2) + 1;
% figure(500); clf; hold on; grid on
% plot(data_x,data_y,'.');


%% estimte the a & b parameters using x & y values
N = numel(data_x);
Z = zeros(N,1);
H = zeros(N,3);

for i = 1:N
   d = log(1/(1 - data_y(i)));
   x = data_x(i);
   x1 = x^2; x2 = 2*x; x3 = 1;
   H(i,:) = [x1 x2 x3];
   Z(i) = d;
end


Xhat = inv(H'*H)*H'*Z;

a = sqrt(Xhat(1));
b = sqrt(Xhat(2));

%%
% clear
inc_x    = 0.01;
data_x   = -4:inc_x:50;
input_x  = data_x';
tru_b = [0.06,-1.2,1]'; % b1,b2,b3
noise_std = 0.0;
data_y  = proposedfuncnoisy(tru_b,input_x,noise_std);

b0   = [0.02,-0.2,0.5]';
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
b_hat_prop = nlinfit(input_x,data_y,@proposedfunc,b0)
% b_hat_log  = nlinfit(input_x,data_y,@logfunc,b0)

sumres_prop = sum(abs(data_y - proposedfunc(b_hat_prop,input_x)))
% sumres_log  = sum(abs(data_y - logfunc(b_hat_log,input_x)))

figure(500); clf; hold on; grid on
plot(data_x,data_y,'.');

%% test on real data
eb_logfunc = 0; % enable/disable log function fitting

data_x   = sz_windows(1:sz_idx);
input_x  = data_x';
data_y  = (min_lin2(1:sz_idx))';
[min_minlin2,idx_min] = min(data_y);
% b0   = [0.02,min_minlin2,1]';
b0   = [0.06,-1.2,1]';
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
b_hat_prop = nlinfit(input_x,data_y,@proposedfunc,b0);
y_propfunc = proposedfunc(b_hat_prop,input_x);
sumres_prop = sum(abs(data_y - y_propfunc));

if eb_logfunc
    b_hat_log  = nlinfit(input_x,data_y,@logfunc,b0);
    y_logfunc  = logfunc(b_hat_log,input_x);
    sumres_log = sum(abs(data_y - y_logfunc));
end

figure(500); clf; hold on; grid on
plot(data_x,data_y,'k.');
plot(data_x,y_propfunc,'-b');
if eb_logfunc
    plot(data_x,y_logfunc,'.r');
    legend({'min. lin.','prop func.','log func.'},'Location','southeast')
else
    legend({'min. lin.','prop func.'},'Location','southeast')
end
title('Minimum Linearity and various function fittings.')
