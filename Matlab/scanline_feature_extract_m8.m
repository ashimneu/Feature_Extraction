function [list] = scanline_feature_extract_m8(scan, scan_pnt_ID)

scan_pnt_ID = scan_pnt_ID(2:end);% first one is the idx of the scan for this run

list(1).pnt_idx = find(scan_pnt_ID == 0);
list(2).pnt_idx = find(scan_pnt_ID == 1);
list(3).pnt_idx = find(scan_pnt_ID == 2);
list(4).pnt_idx = find(scan_pnt_ID == 3);
list(5).pnt_idx = find(scan_pnt_ID == 4);
list(6).pnt_idx = find(scan_pnt_ID == 5);
list(7).pnt_idx = find(scan_pnt_ID == 6);
list(8).pnt_idx = find(scan_pnt_ID == 7);   

list(1).scn_ln = select(scan, list(1).pnt_idx);
list(2).scn_ln = select(scan, list(2).pnt_idx);
list(3).scn_ln = select(scan, list(3).pnt_idx);
list(4).scn_ln = select(scan, list(4).pnt_idx);
list(5).scn_ln = select(scan, list0(5).pnt_idx);
list(6).scn_ln = select(scan, list(6).pnt_idx);
list(7).scn_ln = select(scan, list(7).pnt_idx);
list(8).scn_ln = select(scan, list(8).pnt_idx);

for i = 1:8
    if (~isempty(list(i).pnt_idx))
       [list(i).cor_pc, list(i).ln_pc] = scanline_feature_extract(list(i).scn_ln);
    end
end

end