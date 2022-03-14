function output = hist_point_separation(pc,nbins)
    % nbins - number of histogram plot bins
    N = size(pc,1);
    dist = zeros(N-1,1);
    for j = 1:N-1
        dist(j) = norm(pc(j+1,:)-pc(j,:),2); 
    end
    
    figure(); hold on; grid on;
    h = histogram(dist,nbins);
%     xlim([0 0.3]);

    if nargout == 1
        output = h;
    end
end

