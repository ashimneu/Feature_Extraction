function [laser] = scanline_feature_extract_m8_v2(p,scan, laserID)
    % from scan of multiple lasers, separate x,y,z,range,theta value of points collected for each laser

    scan_pnt_ID = laserID.ID;
    range       = laserID.range;
    theta       = laserID.theta;
    
    % get index of points for each laser
    laser(1).pnt_idx = find(scan_pnt_ID == 0);
    laser(2).pnt_idx = find(scan_pnt_ID == 1);
    laser(3).pnt_idx = find(scan_pnt_ID == 2);
    laser(4).pnt_idx = find(scan_pnt_ID == 3);
    laser(5).pnt_idx = find(scan_pnt_ID == 4);
    laser(6).pnt_idx = find(scan_pnt_ID == 5);
    laser(7).pnt_idx = find(scan_pnt_ID == 6);
    laser(8).pnt_idx = find(scan_pnt_ID == 7);   
    
    % get x,y,z of points in each lasser scanline
    laser(1).scn_ln = select(scan, laser(1).pnt_idx);
    laser(2).scn_ln = select(scan, laser(2).pnt_idx);
    laser(3).scn_ln = select(scan, laser(3).pnt_idx);
    laser(4).scn_ln = select(scan, laser(4).pnt_idx);
    laser(5).scn_ln = select(scan, laser(5).pnt_idx);
    laser(6).scn_ln = select(scan, laser(6).pnt_idx);
    laser(7).scn_ln = select(scan, laser(7).pnt_idx);
    laser(8).scn_ln = select(scan, laser(8).pnt_idx);

    % get laser range of points in each lasser scanline
    laser(1).range = range(laser(1).pnt_idx);
    laser(2).range = range(laser(2).pnt_idx);
    laser(3).range = range(laser(3).pnt_idx);
    laser(4).range = range(laser(4).pnt_idx);
    laser(5).range = range(laser(5).pnt_idx);
    laser(6).range = range(laser(6).pnt_idx);
    laser(7).range = range(laser(7).pnt_idx);
    laser(8).range = range(laser(8).pnt_idx);

    % get laser angle (theta) of points in each lasser scanline
    laser(1).theta = theta(laser(1).pnt_idx);
    laser(2).theta = theta(laser(2).pnt_idx);
    laser(3).theta = theta(laser(3).pnt_idx);
    laser(4).theta = theta(laser(4).pnt_idx);
    laser(5).theta = theta(laser(5).pnt_idx);
    laser(6).theta = theta(laser(6).pnt_idx);
    laser(7).theta = theta(laser(7).pnt_idx);
    laser(8).theta = theta(laser(8).pnt_idx);

    
    for i = p.lasernum
        if (~isempty(laser(i).pnt_idx))
           [laser(i).cor_pc, laser(i).ln_pc, laser(i).cor_pc_index,...
               laser(i).ln_pc_index, laser(i).otherfeatures] = ...
               scanline_feature_extract_v2(p,laser(i).scn_ln);
        end
    end
end