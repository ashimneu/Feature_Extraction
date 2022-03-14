% given: pc,idx_ref_pc,idx_neighbors

num_pnt = size(pc,1);
spc = pc(idx_nei,:);
pnt_ref = pc(idx_ref_pc,:);

figure(9595); clf; hold on; grid on
plot3(pc(:,1),pc(:,2),pc(:,3),'b.')
plot3(pnt_ref(1),pnt_ref(2),pnt_ref(3),'rp')
% plot3(spc(:,1),spc(:,2),spc(:,3),'c+')
[sorted_idx_nei,sorting_idx] = sort_index(idx_nei);
aspc = spc(sorting_idx,:);
plot3(aspc(:,1),aspc(:,2),aspc(:,3),'c')

% i_pnt_ref_sw = find(idx_nei == idx_ref_pc);
%%
[asc_idx_nei,asc_idx] = sort(idx_nei,'ascend');
idx_n_shifted = shift_index_v2(asc_idx_nei);

[~,inv_asc_idx] = sort(asc_idx,'ascend');

%%
