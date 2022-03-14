function sz_window = checkWindowSize(sz_window,num_pnts,sz_min_max)
    sz_window_min = sz_min_max(1);
    sz_window_max = sz_min_max(2);

    if (sz_window > num_pnts)
        % disp('Sliding window size is larger than the number of points.')
        % disp('Setting window size to max. number of points.')
        sz_window = sz_window_max;
    elseif (sz_window < sz_window_min)
        % disp('Sliding window size is less than threshold of 5.')
        % disp('Setting window size to 5.')
        sz_window = sz_window_min;
    end
end

