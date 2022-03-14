function output = hist_rangetheta_separation(range,theta,nbins)
    % nbins - number of histogram plot bins

    rdiff = abs(range(2:end) - range(1:end-1));
    tdiff = theta(2:end) - theta(1:end-1);
    
    figure();
    subplot(211); hold on; grid on;
    h1 = histogram(rdiff,nbins);
    title('range')
    subplot(212); hold on; grid on;
    h2 = histogram(tdiff,nbins);
    title('range')

    if nargout == 1
        output.h1 = h1;
        output.h2 = h2;
    end
end
