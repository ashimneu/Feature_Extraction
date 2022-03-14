function append_DataTipTemplate(tipname,gs,cluster_id,idx_cor_L)
    [uniq_id,~,~] = unique(cluster_id);
    for clust_id = 1:max(uniq_id)
        % index of id-th cluster points in pnt_cor_L (pnt_cor_L is set 
        % of corner points extracted using linearity check).
        idxa = find(cluster_id == clust_id); 
        % index of id-th cluster points in pc_xyz (pc_xyz is pc of single laser)
        idxb = idx_cor_L(idxa);  
        gs(clust_id).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow(tipname,idxb);
    end
end

