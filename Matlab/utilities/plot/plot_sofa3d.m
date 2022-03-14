function plot_sofa3d(sz_idx,sz_windows,sofa_3d)
sz_windows = sz_windows(1:sz_idx);
for i = 1:sz_idx
    window = sz_windows(i);
    SOFA = sofa_3d{i};
    num = numel(SOFA);
    x = window.*ones(num,1);
    y = 1:num;
    plot3(x,y,SOFA,'.')
end

end

