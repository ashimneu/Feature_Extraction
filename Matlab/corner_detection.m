% Note: Please run this script (line_detection.m) only after
% running process_scanline_new_features.m

% Variables required from process_scanline_new_features.m are
% pc_xyz, cor_pc, otherfeatures

% Acronyms:
% ln   - line; L-linearity; P-Planarity; S-Sphericity; O-Omnivariance;
% Esum - EigenValue Sum; Entropy-EigenEntropy; Sm-Smoothness; dp-DotProduct
% SV   - surface variation; A-Anisotropy
% FOFA - First Order, First Axis Moment;  FOSA - First Order, Second Axis Moment
% SOFA - Second Order, First Axis Moment; SOSA - Second Order, Second Axis Moment
% pnt_clust - points in current cluster

clc;

%parameters
p.s = rng; % save current random number generator seed
p.sw_def = 3; % 2,3 are knnSearch v2 & distance based window, respectively 
p.threshold_L = 0.85;
p.threshold_dp = 0.95;

% observe histogram to gauge max_dist value between points in a cluster
diff = cor_pc(2:end,:) - cor_pc(1:end-1,:);
diff_norm = vecnorm(diff,2,2);
% enable this histogram to gauge epsilon from pc data
% figure(201); hold on; grid on;
% h = histogram(diff_norm,100); 

% retrieve the features computed in process_scanline_new_features.m
[L,P,S,O,A,Esum,SV,Entropy,Sm,FOFA,FOSA,SOFA,SOSA,dp,angle] = deal(otherfeatures{:});

% Based on Linearity check, get corner points & index in the scanline 
[idx_cor_L,~,pnt_cor_L,~] = find_index(L,p.threshold_L,pc_xyz); 

% DBSCAN to extract corner cluster & filter out outliers
minpts = 5;   % min points for dbscan
max_dist  = 1; % max distance between points in a cluster
cluster_id = dbscan(pnt_cor_L,max_dist,minpts); % get id of unique cluster points

% Plot colored clusters of the extracted corner points
plot_clusters(202,pnt_cor_L,cluster_id,idx_cor_L); % fignum:202

[uniq_id,ia,ic] = unique(cluster_id);
uniq_id(1) = []; % discard id represnting outliers (i.e. id = -1 is discarded).
idxa = [];
idxb = [];

sz_windows = 5:2:38; % 38, 30
radiuses = .4:.2:.8;
min_L2 = nan(numel(sz_windows));
sofa_3d  = cell(numel(sz_windows),1);

max_FOSA2_idx = nan(numel(radiuses),1);
min_FOSA2_idx = nan(numel(radiuses),1);
max_SOSA2_idx = nan(numel(radiuses),1);
min_SOSA2_idx = nan(numel(radiuses),1);
min_dp2_idx   = nan(numel(radiuses),1);

for sz_idx = 1:numel(radiuses)
% for sz_idx = 1:numel(sz_windows)
for clust_id = 7%1:max(uniq_id)
    idxa = find(cluster_id == clust_id); % index of id-th cluster points in pnt_cor_L (pnt_cor_L is set of corner points extracted using linearity check).
    pnt_clust = pnt_cor_L(idxa,:);
    [medoid,~] = mfames(pnt_clust'); % medoid
    idxb = idx_cor_L(idxa);  % index of id-th cluster points in pc_xyz (pc_xyz is pc of single laser)
    num_pnt = numel(idxa);
    len_clust = computelength(pnt_clust); % length of cluster

    sz_window = 11; % sz_windows(sz_idx);
    p.radius = radiuses(sz_idx);

    fprintf('radius = %2.2f. \n',p.radius)
    % fprintf('sliding window size = %2.0f. \n',sz_window)

    % recompute Geometric Features on the idth cluster
    %     L2 = compute_linearity_using_sw(pnt_clust,sz_window);    
    %     O2 = compute_omnivar_using_sw(pnt_clust,sz_window);
    %     [FOFA2,FOSA2,SOFA2,SOSA2] = compute_moment_using_sw(pnt_clust,sz_window);
    features2 = compute_feature_using_sw_v2(p,pnt_clust,sz_window);
    [L2,P2,S2,O2,A2,Esum2,SV2,Entropy2,Sm2,FOFA2,FOSA2,SOFA2,SOSA2,dp2,len_pr2,len_sc2] = deal(features2{:});
    sofa_3d{sz_idx} = SOFA2;
    min_L2(sz_idx) = min(L2); % minimum linearity when using jth sw_size 
    
    % Find corner point using EigenSum
    [maxESum, idx_maxESum] = max(Esum(idxb));
    pnt_maxESum = pc_xyz(idxb(idx_maxESum),:); % point corresponding to maxEigSum
    
    % Find corner point using dot product of fitted lines
    [mindp, idx_mindp] = min(dp2);
    pnt_mindp = pc_xyz(idxb(idx_mindp),:); % point corresponding to maxEigSum
   
    % Plot features computed on entire scanline
    gr = 3; gc = 2; % subplot grid rows & cols count
    figure(301); clf;
    subplot(gr,gc,1); hold on; grid on
    s1 = scatter3(pnt_clust(:,1),pnt_clust(:,2),pnt_clust(:,3),'.');
    s1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb);
    % plot3(pnt_maxESum(1),pnt_maxESum(2),pnt_maxESum(3),'rp')
    xlabel('x'); ylabel('y'); zlabel('z');
    title({['Cluster ',num2str(clust_id),' point cloud ']});
    % axis equal

    subplot(gr,gc,2); hold on; grid on
    pnt_clust_trans = pnt_clust - medoid';
    s1 = scatter3(pnt_clust_trans(:,1),pnt_clust_trans(:,2),pnt_clust_trans(:,3),'.');
    s1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idxb);
    % plot3(pnt_maxESum(1),pnt_maxESum(2),pnt_maxESum(3),'rp')
    xlabel('x'); ylabel('y'); zlabel('z');
    title({['Translated Cluster ',num2str(clust_id),' point cloud ']});

    subplot(gr,gc,3); hold on; grid on
    plot(idxb,Esum(idxb),'.');
    ylabel('EigSum')
    
    subplot(gr,gc,4); hold on; grid on
    plot(idxb,L(idxb),'.');
    ylabel('Linearity')
    
    subplot(gr,gc,5); hold on; grid on
    plot(idxb,dp(idxb),'.');
    ylabel('Dot Product')
    
    subplot(gr,gc,6); hold on; grid on
    plot(idxb,angle(idxb),'.');
    ylabel('Angle')

    %     subplot(gr,gc,4); hold on; grid on
    %     plot(idxb,O(idxb),'.');
    %     ylabel('Omnivariance')

    % subplot(gr,gc,5); hold on; grid on
    % plot(idxb,SV(idxb),'.');
    % ylabel('Surf. Var.')

    % subplot(gr,gc,6); hold on; grid on
    % plot(idxb,Sm(idxb),'.');
    % ylabel('Smoothness')

    % subplot(gr,gc,7); hold on; grid on
    % plot(idxb,FOSA(idxb),'.');
    % ylabel('FOSA')
    % xlabel('index')
    % 
    % subplot(gr,gc,8); hold on; grid on
    % plot(idxb,SOSA(idxb),'.');
    % ylabel('SOSA')
    % xlabel('index')

    sgtitle({['Features computed on laser #',num2str(lzid),' pointcloud.'], ...
        ['cluster #',num2str(clust_id),', point count=',num2str(num_pnt),...
        ', len=',num2str(round(len_clust,2)),'m']});

    
    % Plot the features computed on idth cluster
    gr = 3; gc = 2; % subplot grid rows & cols count
    fi = 0;
    % Plot 1st & 2nd moments
    figure(302); clf; 
    
%     fi = fi+1;
%     subplot(gr,gc,fi);hold on; grid on
%     plot(idxb,FOFA2,'.')
%     ylabel('FOFA')
%     title({['cluster ',num2str(clust_id)]});
    
    subplot(gr,gc,1);hold on; grid on
    plot(idxb,FOSA2,'.')
    ylabel('FOSA')
    [~,i_maxFOSA2] = max(FOSA2);
    [~,i_minFOSA2] = min(FOSA2);
    title({['max:',num2str(idxb(i_maxFOSA2)),' min:',num2str(idxb(i_minFOSA2))]});
    max_FOSA2_idx(sz_idx) = idxb(i_maxFOSA2);
    min_FOSA2_idx(sz_idx) = idxb(i_minFOSA2);

%     subplot(gr,gc,3);hold on; grid on
%     plot(idxb,SOFA2,'.')
%     ylabel('SOFA')
    
    subplot(gr,gc,2);hold on; grid on
    plot(idxb,SOSA2,'.')
    ylabel('SOSA');
    [~,i_maxSOSA2] = max(SOSA2);
    [~,i_minSOSA2] = min(SOSA2);
    title({['max:',num2str(idxb(i_maxSOSA2)),' min:',num2str(idxb(i_minSOSA2))]});
    max_SOSA2_idx(sz_idx) = idxb(i_maxSOSA2);
    min_SOSA2_idx(sz_idx) = idxb(i_minSOSA2);

    subplot(gr,gc,3); hold on; grid on
    plot(idxb,L2,'.');
    ylabel('Linearity')
    xlabel('index')

    subplot(gr,gc,4); hold on; grid on
    plot(idxb,dp2,'.');
    ylabel('dot product')
    xlabel('index')
    [~,i_mindp2] = min(dp2);
    title({['min:',num2str(idxb(i_mindp2))]});
    min_dp2_idx(sz_idx) = idxb(i_mindp2);

    subplot(gr,gc,5); hold on; grid on
    plot(idxb,Esum2,'.');
    ylabel('EigenSum')
    xlabel('window size')

    subplot(gr,gc,6); hold on; grid on
    plot(idxb,len_pr2,'r-+');
    plot(idxb,len_sc2,'c-+');
    ylabel('fitted line length')
    xlabel('index')
    legend('pr','sc','location','south')
    

    if p.sw_def == 1 || p.sw_def == 2
        sgtitle({['Features computed on cluster #',num2str(clust_id),' only.'], ...
            ['window size=',num2str(sz_window), ', point count=',num2str(num_pnt)]});
    elseif p.sw_def == 3
        sgtitle({['Features computed on cluster #',num2str(clust_id),' only.'], ...
                ['radius=',num2str(p.radius),'m, point count=',num2str(num_pnt)]});
    end
    
    %     if sz_idx == 1
    %         eigsum_cluster = Esum(idxb);
    %         threshold = 0.02;
    %         
    %         [fit_flag,cost_1storderfit,cost_2ndorderfit] = lin_quad_fit_costcompare(idxb,eigsum_cluster,threshold);
    %         % linear fit: flag=1; nonlinear fit: flag = 2;
    %         fprintf('Ploynomial fitting on Eigen value Sum \n')
    %         fprintf('plot (see fig. 301) is performed. \n')
    %         fprintf('1st order fit cost = %2.4f. \n',cost_1storderfit);
    %         fprintf('2nd order fit cost = %2.4f. \n',cost_2ndorderfit);
    %         if fit_flag == 1
    %             fprintf('Points in cluster%1.0f are linear. No corner detected.\n',clust_id);
    %         elseif fit_flag == 2
    %             fprintf('Points in cluster%1.0f are non-linear.\n',clust_id);
    %             fprintf('Index of corner point= %6.0f \n',idxb(idx_maxESum));
    %         end
    %     end
        
    
    
    % compute angle between fitted lines using FOSA based corner detection
    modeFOSA2 = mode(max_FOSA2_idx); % idx_maxESum;
    idx_mode_FOSA2 = find(idxb == modeFOSA2); 
    [~,theta_FOSA,~,~] = dotprod_fitlines_v1(p,pnt_clust,idx_mode_FOSA2);

    % compute angle between lines using dotproduct based corner detection
    modedp2 = mode(min_dp2_idx); % idx_maxESum;
    idx_mode_dp2 = find(idxb == modedp2); 
    [~,theta_dp,~,~] = dotprod_fitlines_v1(p,pnt_clust,idx_mode_dp2);
    
    %     pr_pnts = pnt_clust(1:idx_mode,:);   % points preceeding the one having maxEigSum (including maxEigSum point)
    %     sc_pnts = pnt_clust(idx_mode:end,:); % points succeeding the one having maxEigSum (including maxEigSum point)
    %         
    %     % use svd to get vector of line fit on preeceding & succedding points
    %     [U_pr,sv_pr,~] = svd(pr_pnts');
    %     [U_sc,sv_sc,~] = svd(sc_pnts');
    %     fprintf('\nSingular values of \n')
    % %     fprintf('preeceding points: %3.3f, %3.3f, %3.3f \n',[sv_pr(1,sv_pr(2,2),sv_pr(3,3)]);
    % %     fprintf('succeeding points: %3.3f, %3.3f, %3.3f \n',diag(sv_sc));   
    % 
    %     % compute incident angle between the two lines
    %     cosangle = @(u,v) dot(u,v)/(norm(u)*norm(v)); % cosine of incident angle
    %     incangle = @(cosangle) acosd(cosangle); % incident in degrees
    % 
    %     cosalpha = cosangle(U_pr(:,1),U_sc(:,1));
    %     alpha = incangle(cosalpha);
    %     fprintf('incident angle: %3.3f deg. \n',alpha);
end
end