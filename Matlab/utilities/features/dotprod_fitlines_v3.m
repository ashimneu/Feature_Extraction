function [dotprod,alpha,len_pr,len_sc] = dotprod_fitlines_v3(p,pc,pnt_ref_idx)
    % for a given point cloud and a reference point in the point cloud,
    % point cloud is divided into two halves. First half are the preceeding
    % points from reference point and second half are succeeding points.
    % Stright line is fitted on both halves and unit vector along both 
    % the fitted lines is computed. Dot product of the unit vector is
    % computed here as a feature to detect corner in the point cloud.
    % The input pc is not feasible for the number of points is less than
    % minimum number of points required for line fitting.
    % INPUT: 
    % pnt_ref_idx - index in pc array of a reference point for which
    %               the feature is computed.
    % pc - point cloud
    
    num_pnt = size(pc,1);   
  
    % Check if reference point (the one for which dot product is computed) 
    % is on either end of pc, there will be insufficient preceeding or
    % succeeding points required for line fitting.
    if  pnt_ref_idx == 1 || pnt_ref_idx == 2 || pnt_ref_idx == 3 || pnt_ref_idx == (num_pnt-2) || pnt_ref_idx == (num_pnt-1) || pnt_ref_idx == num_pnt
    % if  pnt_ref_idx == 1 || pnt_ref_idx == 2 || pnt_ref_idx == (num_pnt-1) || pnt_ref_idx == num_pnt
    % if  pnt_ref_idx == 1 || pnt_ref_idx == num_pnt
        idx_is_invalid = true;
    else
        idx_is_invalid = false;
    end
    
    % unused variables
    theta   = nan;
    len_pr  = nan;
    len_sc  = nan;

    if idx_is_invalid
        dotprod = 3;
        alpha   = -3;        
    else
        % reference point = medoid or centroid or center point in pc array
        pnt_ref = pc(pnt_ref_idx,:); 
        
        % translate the pc such that reference point is at the origin
        translated_pc = pc - pnt_ref;
        pr_pnts = translated_pc(1:pnt_ref_idx,:);   % <S1>, preceeding points
        sc_pnts = translated_pc(pnt_ref_idx:end,:); % <S2>, succeeding points
        
        % svd over preceeding points (step 2a)
        [U_pr,S_pr,~] = svd(pr_pnts');
        sv_pr1 = S_pr(1,1); sv_pr2 = S_pr(2,2);
        U_pr1  = U_pr(:,1); U_pr2  = U_pr(:,2);

        % svd over succeeding points (step 2b)
        [U_sc,S_sc,~] = svd(sc_pnts');
        sv_sc1 = S_sc(1,1); sv_sc2 = S_sc(2,2);
        U_sc1 = U_sc(:,1); U_sc2 = U_sc(:,2);
        
        % (step 3) check if preceeding & succeeding points are linear points
        sv_threshold = 2; % singular value threshold to satisfy single dimensionality of the points
        pr_ratio = (sv_pr1 - sv_pr2 )/sv_pr2;
        sc_ratio = (sv_sc1 - sv_sc2 )/sv_sc2;
        pr_is_linear = pr_ratio > sv_threshold; % threshold check
        sc_is_linear = sc_ratio > sv_threshold; % threshold check

        % PCA over preceeding points (step 2a)
        [V_pr,E_pr,~] = compute_pca(pr_pnts');
        E_pr1 = E_pr(1,1); E_pr2 = E_pr(2,2);
        V_pr1 = V_pr(:,1); V_pr2 = V_pr(:,2);

        % PCA over succeeding points (step 2b)
        [V_sc,E_sc,~] = compute_pca(sc_pnts');
        E_sc1 = E_sc(1,1); E_sc2 = E_sc(2,2);
        V_sc1 = V_sc(:,1); V_sc2 = V_sc(:,2);

        % (step 3) check if preceeding & succeeding points are linear points
        E_threshold = 10; % singular value threshold to satisfy single dimensionality of the points
        pr_ratio2 = E_pr1/E_pr2;
        sc_ratio2 = E_sc1/E_sc2;
        pr_is_linear2 = pr_ratio2 > E_threshold; % threshold check
        sc_is_linear2 = sc_ratio2 > E_threshold; % threshold check
        
        % dotprod = abs(dot(U_pr1,U_sc1)); % dot product    
        % (step 4) compute dot product if preceeding & succedding points are linear
        if pr_is_linear2 && sc_is_linear2
            dotprod2 = abs(dot(U_pr1,U_sc1)); % dp using svd only
            dotprod = abs(dot(V_pr1,V_sc1)); % dot product after PCA
            cos_alpha = dot(V_pr1,V_sc1)/(norm(V_pr1)*norm(V_sc1));            
            alpha = 180 - acosd(max(-1,min(cos_alpha,1))); % angle between lines, min & max because -1 < cos(theta) < 1 for all real theta.
        else
            dotprod = 2;
            alpha   = -2;
        end

        % indices for debugging 
        debug_indices = []; %indices; <--- INSERT INDEX IF DEBUGGING ELSE KEEP ARRAY EMPTY ------

        do_debug = sum(p.idx_debug == debug_indices);
        if  do_debug
            figure(); hold on; grid on
            plt_pr = plot3(pr_pnts(:,1),pr_pnts(:,2),pr_pnts(:,3),'.r');
            plt_sc = plot3(sc_pnts(:,1),sc_pnts(:,2),sc_pnts(:,3),'.b');
            plt_pr.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',p.idxa_sw(1:pnt_ref_idx));
            plt_sc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',p.idxa_sw(pnt_ref_idx:end));
            plot3(0,0,0,'p')
            axis equal
            legend('pr','sc','location','best')
            xlabel('x'); ylabel('y'); zlabel('z')
            title({['for index ',num2str(p.idx_debug)]})
            
            fprintf('\nIndex: %4.0f \n',p.idx_debug)
            fprintf('PCA: \n')
            fprintf('E1_pr = %1.3f, E2_pr = %1.3f, rat2_pr = %3.2f \n',[E_pr1,E_pr2,pr_ratio2])
            fprintf('E1_sc = %1.3f, E2_sc = %1.3f, rat2_sc = %3.2f \n',[E_sc1,E_sc2,sc_ratio2])
            fprintf('dotprod = %1.3f \n',dotprod)
            fprintf('angle = %1.3f \n',alpha)
            [pr_is_linear2,sc_is_linear2]
            direction2 = [V_pr1';V_sc1']
            direction2 = [V_pr1';V_sc1']; % <-----PUT LINE BREAK HERE FOR DEBUG-------
        end
        
        % svd over translated pc
        [U,E,~] = svd(translated_pc');
        u1 = U(:,1); u2 = U(:,2); u3 = U(:,3); % eigen vectors
        e1 = E(1,1); e2 = E(2,2); e3 = E(3,3); % eigen values
        
        % %%%%% Please ignore the codes after this!!
        % % u-direction, v-direction
        % pr_u1 = u1'*pr_pnts'; % u component of pr
        % pr_u2 = u2'*pr_pnts'; % v component of pr
        % sc_u1 = u1'*sc_pnts'; % u component of sc
        % sc_u2 = u2'*sc_pnts'; % v component of sc    
        % 
        % % MSAC (robust variant of RANSAC)
        % sampleSize  = 3;   % number of points to sample per trial
        % maxDistance = 0.5; % max allowable distance for inliers
        % 
        % fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
        % % distance evaluation function
        % evalLineFcn = @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
        % 
        % rng(p.s); % use fixed rng seed so that ransac's result doesn't vary. 
        %           % Ransac computation uses randomly selected intial points so
        %           % its results aren't always same for idential inputs.
        %     
        % points_pr = [pr_u1', pr_u2'];
        % points_sc = [sc_u1', sc_u2'];
        % 
        %
        % algo_linfit = 2; % line fitting algorithm: 1-ransac, 2-ols
        % switch algo_linfit
        %     case 1 % RANSAC
        %         [~, inlierIdx_pr] = ransac(points_pr,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
        %         model_pr = polyfit(pr_u1(inlierIdx_pr)',pr_u2(inlierIdx_pr)',1);
        %         yhat_pr = polyval(model_pr,pr_u1);
        %         [~, inlierIdx_sc] = ransac(points_sc,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
        %         model_sc = polyfit(sc_u1(inlierIdx_sc)',sc_u2(inlierIdx_sc)',1);
        %         yhat_sc = polyval(model_sc,sc_u1);
        %     case 2 % Ordinary Least Squares
        %         model_pr = polyfit(pr_u1',pr_u2',1);
        %         yhat_pr = polyval(model_pr,pr_u1);
        %         model_sc = polyfit(sc_u1',sc_u2',1);
        %         yhat_sc = polyval(model_sc,sc_u1);            
        % end
        % 
        % % compute length of both fitted lines
        % len_pr = norm([pr_u1(1); yhat_pr(1)] - [pr_u1(end); yhat_pr(end)]);
        % len_sc = norm([sc_u1(1); yhat_sc(1)] - [sc_u1(end); yhat_sc(end)]);
        % 
        % % compute angle between fitted lines using normal to the lines
        % n_pr = [1 -model_pr(1)]; % normal vector = coefficients of x & y of line
        % n_sc = [1 -model_sc(1)]; % normal vector = coefficients of x & y of line
        % costheta = dot(n_pr,n_sc)/(norm(n_pr)*norm(n_sc));
        % theta = 180 - acosd(costheta); % angle between the fitted lines
        % 
        % vec_pr = [pr_u1(1), yhat_pr(1)] - [pr_u1(end), yhat_pr(end)];
        % vec_sc = [sc_u1(1), yhat_sc(1)] - [sc_u1(end), yhat_sc(end)];
        % vec_pr = vec_pr/norm(vec_pr);
        % vec_sc = vec_sc/norm(vec_sc);
        % % dotprod2 = abs(dot(vec_pr,vec_sc)); % dot product
        % 
        % DEBUG = false;
        % if DEBUG
        %     
        %     figure(501); clf; 
        %     subplot(211); hold on; grid on        
        %     
        %     % plot eigen vectors
        %     plot3([0 u1(1)],[0 u1(2)],[0 u1(3)],'-m')
        %     plot3([0 u2(1)],[0 u2(2)],[0 u2(3)],'-g')
        %     % plot preceeding & succedding points
        %     plt_pr = plot3(pr_pnts(:,1),pr_pnts(:,2),pr_pnts(:,3),'.r');
        %     plt_sc = plot3(sc_pnts(:,1),sc_pnts(:,2),sc_pnts(:,3),'.c');
        %     plt_pr.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',p.idxa_sw(1:pnt_ref_idx));
        %     plt_sc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',p.idxa_sw(pnt_ref_idx:end));
        %     legend('pr','sc','u1','u2')
        %     xlabel('x'); ylabel('y'); zlabel('z')
        %     % axis equal
        %     
        %     subplot(212); hold on; grid on
        %     plot(pr_u1,pr_u2,'.r')
        %     plot(sc_u1,sc_u2,'.c')
        %     xlabel('u1'); ylabel('u2')
        %     
        %     subplot(212); hold on; grid on
        %     plot(pr_u1,yhat_pr,'-r')
        %     plot(sc_u1,yhat_sc,'-c')
        %     title('projected along u1 & u2')
        %     legend('pr','sc','linefit on pr','linefit on sc')
        % end
    end
end

