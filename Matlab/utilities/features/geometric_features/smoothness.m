function Sm = smoothness(pc,pnt_ref_idx)
    % pc_scanline - (Nx3) points of a single laser scan
    % pnt_ref_ids - index of the point for which smoothness is computed.

    num_pnts = size(pc,1);
    pnt_ref  = pc(pnt_ref_idx,:);    
    pc(pnt_ref_idx,:) = []; % note: reference point is excluded
    diff     = pc - pnt_ref;
    diffsum  = sum(diff,1);
    diffsum_norm = norm(diffsum);
    ref_norm = norm(pnt_ref);
    Sm = diffsum_norm / (num_pnts * ref_norm);
end
