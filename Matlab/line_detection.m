% Note: Please run this script (line_detection.m) only after
% running process_scanline_new_features.m

% Variables required from process_scanline_new_features.m are
% pc_xyz, ln_pc, otherfeatures

% Acronyms:
% ln   - line
% SV   - surface variation
% FOFA - First Order, First Axis Moment
% FOSA - First Order, Second Axis Moment
% SOFA - Second Order, First Axis Moment
% SOSA - Second Order, Second Axis Moment

clc;

% observe histogram to gauge max_dist value between points in a cluster
diff = cor_pc(2:end,:) - cor_pc(1:end-1,:);
diff_norm = vecnorm(diff,2,2);
% enable this histogram to gauge epsilon from pc data
% figure(201); hold on; grid on;
% h = histogram(diff_norm,100); 

[L,P,S,O,A,Esum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA] = deal(otherfeatures{:});

% Based on Linearity check, get identify stright line points & index in the scanline 
% pnt_cor_L & idx_ln_L are points detected as corner & their indices, respectively.
[~,idx_ln_L,~,pnt_ln_L] = find_index(L,p.threshold_L,pc_xyz); 
idx_endpoint = max(idx_ln_L); % index of last data in the scanline pointCloud.

% DBSCAN to extract corner cluster & filter out outliers
minpts     = 10;  % min points for dbscan
max_dist   = 0.3; % 0.5; % max distance between points in a cluster
cluster_id = dbscan(pnt_ln_L,max_dist,minpts); % get id of unique cluster points

% Plot colored clusters of the extracted line points
plot_clusters(202,pnt_ln_L,cluster_id); % fignum:202

[uniq_id,ia,ic] = unique(cluster_id);
uniq_id(1) = []; % discard id represnting outliers (i.e. id = -1 is discarded).
idxa = [];
idxb = [];

for clust_id = 1 % 1:max(uniq_id)

    idxa = find(cluster_id == clust_id); % index of id-th cluster points in the array of corner points (i.e. pnt_cor_L) 
    pnt_clust = pnt_ln_L(idxa,:); % points in the cluster
    idxb = idx_ln_L(idxa); % index of id-th cluster points in the entire scanline (i.e. pc_xyz)
    
    % check if cluster contains leading & trailing ends of the scanline.
    % if yes, to avoid index dropping from end to 1 during feature computation
    % reassign index of all points starting from index = 1
    if  contains_ends(idxb,idx_endpoint)
        idxb2 = shift_index(idxb);
        [idxb2,idx_sort] = sort(idxb2,'ascend');
        pnt_clust2 = pnt_clust(idx_sort,:);
        idxb = idxb(idx_sort);
    else
        idxb2 = idxb;
        pnt_clust2 = pnt_clust;
    end

    num_pnt    = size(pnt_clust,1); % number of points in the cluster
    sz_windows = 5:2:odd_int(num_pnt); % 38, 30 % list of window size to evalute
    
    % Allocate variables for features storage
    min_lin2 = nan(numel(sz_windows));
    sofa_3d  = cell(numel(sz_windows),1);

for sz_idx = 1:numel(sz_windows)
    sz_window = sz_windows(sz_idx); % sliding window size
    fprintf('sliding window size = %2.0f. \n',sz_window)
    
    % compute geometric features on idth cluster
    linearity2 = compute_linearity_using_sw(pnt_clust2,sz_window);
    min_lin2(sz_idx) = min(linearity2); % minimum linearity when using jth sw_size 
    omnivariance2 = compute_omnivar_using_sw(pnt_clust2,sz_window);
    [FOFA2,FOSA2,SOFA2,SOSA2] = compute_moment_using_sw(pnt_clust2,sz_window);
    sofa_3d{sz_idx} = SOFA2;
    [maxEigSum, idx_maxEigSum] = max(Esum(idxb));

    gr = 4; gc = 3; % subplot grid rows & cols count
    figure(301);clf;
    subplot(gr,gc,1); hold on; grid on
    s1 = scatter3(pnt_clust2(:,1),pnt_clust2(:,2),pnt_clust2(:,3),'.');
    s1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb);
    xlabel('x'); ylabel('y'); zlabel('z');
    title({['Cluster ',num2str(clust_id),' point cloud ']});
    axis equal

    subplot(gr,gc,2); hold on; grid on
    plot(idxb2,L(idxb),'.');
    ylabel('Linearity')

    subplot(gr,gc,3); hold on; grid on
    plot(idxb2,Esum(idxb),'.');
    ylabel('EigSum')

    subplot(gr,gc,4); hold on; grid on
    plot(idxb2,O(idxb),'.');
    ylabel('Omnivariance')

    subplot(gr,gc,5); hold on; grid on
    plot(idxb2,SV(idxb),'.');
    ylabel('Surf. Var.')

    subplot(gr,gc,6); hold on; grid on
    plot(idxb2,Sm(idxb),'.');
    ylabel('Smoothness')

    subplot(gr,gc,7); hold on; grid on
    plot(idxb2,FOFA(idxb),'.');
    ylabel('FOFA')
    xlabel('index')

    subplot(gr,gc,8); hold on; grid on
    plot(idxb2,FOSA(idxb),'.');
    ylabel('FOSA')
    xlabel('index')

    subplot(gr,gc,9); hold on; grid on
    plot(idxb2,SOFA(idxb),'.');
    ylabel('SOFA')
    xlabel('index')

    subplot(gr,gc,10); hold on; grid on
    plot(idxb2,SOSA(idxb),'.');
    ylabel('SOSA')
    xlabel('index')

    sgtitle({['Features computed on laser #',num2str(lzid),' pointcloud.'], ...
        ['cluster #',num2str(clust_id),', point count=',num2str(num_pnt)]});

    
    gr = 4; gc = 2; % subplot grid rows & cols count
    % Plot 1st & 2nd moments
    figure(302); clf; 
    subplot(gr,gc,1);hold on; grid on
    p1 = plot(idxb2,FOFA2,'.');
    ylabel('FOFA')
    title({['cluster ',num2str(clust_id)]});    
    
    subplot(gr,gc,2);hold on; grid on
    p2 = plot(idxb2,FOSA2,'.');
    ylabel('FOSA')
    
    subplot(gr,gc,3);hold on; grid on
    p3 = plot(idxb2,SOFA2,'.');
    ylabel('SOFA')
    
    subplot(gr,gc,4);hold on; grid on
    p4 = plot(idxb2,SOSA2,'.');
    ylabel('SOSA');
%     axis equal

    subplot(gr,gc,5); hold on; grid on
    p5 = plot(idxb2,linearity2,'.');
    ylabel('Linearity')
    xlabel('index')

    subplot(gr,gc,6); hold on; grid on
    p6 = plot(idxb2,omnivariance2,'.');
    ylabel('Omnivariance')
    xlabel('index')
    
    subplot(gr,gc,7); hold on; grid on
    plot(sz_windows(1:sz_idx),min_lin2(1:sz_idx),'.');
    ylabel('min linearity')
    xlabel('window size')
    
    subplot(gr,gc,8); hold on; grid on
    plot_sofa3d(sz_idx,sz_windows,sofa_3d);
    zlabel('SOFA')
    ylabel('index')
    xlabel('window size')
    title('SOFA vs window size')

    sgtitle({['Features computed on cluster #',num2str(clust_id),' only.'], ...
        ['window size=',num2str(sz_window), ', point count=',num2str(num_pnt)]});
    
    p1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb(idx_sort));
    p2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb(idx_sort));
    p3.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb(idx_sort));
    p4.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb(idx_sort));
    p5.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb(idx_sort));
    p6.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb(idx_sort));

    if sz_idx == 1
    eigsum_cluster = Esum(idxb);
    threshold = 0.05;
    [fit_flag,cost_1storderfit,cost_2ndorderfit] = lin_quad_fit_costcompare(idxb,eigsum_cluster,threshold);
    % linear fit: flag=1; nonlinear fit: flag = 2;
    fprintf('Ploynomial fitting on Eigen value Sum \n')
    fprintf('plot (see fig. 301) is performed. \n')
    fprintf('1st order fit cost = %2.4f. \n',cost_1storderfit);
    fprintf('2nd order fit cost = %2.4f. \n',cost_2ndorderfit);
    fprintf('fit_flag %1.0f \n',fit_flag);
%     if fit_flag == 1
%         fprintf('Points in cluster%1.0f are linear. No corner detected.\n',clust_id);
%     elseif fit_flag == 2
%         fprintf('Points in cluster%1.0f are non-linear.\n',clust_id);
%         fprintf('Index of corner point= %6.0f \n',idxb(idx_maxEigSum));
%     end
    end
end
end