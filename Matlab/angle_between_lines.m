function theta = angle_between_lines(pc,idx)
    % pnt_clust - cluster points
    % pr_pnt - preceeding points in the cluster
    % sc_pnt - succeeding points in the cluster
    
    pr_pnts = pc(1:idx,:);
    sc_pnts = pc(idx:end,:);
    
    % [med,~ ] = mfames(pnt_clust'); % medoid
    % cent = mean(pc_xyz,1)'; % centroid
    % const = med;
    const = [0; 0; 0];
    pnt1 = pr_pnts - const';
    pnt2 = sc_pnts - const';
    
    % % plot preceeding & succedding points
    % figure(501); clf; 
    % subplot(211); hold on; grid on
    % plot3(pnt1(:,1),pnt1(:,2),pnt1(:,3),'.r')
    % plot3(pnt2(:,1),pnt2(:,2),pnt2(:,3),'.c')
    % axis equal
    
    % pr_pnts; sc_pnts
    [U,E,~] = svd(pc');
    u1 = -U(:,1); u2 = -U(:,2); u3 = U(:,3); % eigen vectors
    e1 = E(1,1); e2 = E(2,2); e3 = E(3,3); % eigen values
    
    % % plot eigen vectors
    % plot3([0 u(1)],[0 u(2)],[0 u(3)],'-k')
    % plot3([0 v(1)],[0 v(2)],[0 v(3)],'-m')
    % xlabel('x'); ylabel('y'); zlabel('z')
    % legend('pr','sc','u1','u2')
    % % axis equal
    
    % % u-direction, v-direction
    % pr_u1 = u1'*pr_pnts'; % u1 component of pr
    % pr_u2 = u2'*pr_pnts'; % u2 component of pr
    % sc_u1 = u1'*sc_pnts'; % u1 component of sc
    % sc_u2 = u2'*sc_pnts'; % u2 component of sc
    
    % u-direction, v-direction
    pr_u1 = u1'*pnt1'; % u component of pr
    pr_u2 = u2'*pnt1'; % v component of pr
    sc_u1 = u1'*pnt2'; % u component of sc
    sc_u2 = u2'*pnt2'; % v component of sc
    
    
    % subplot(212); hold on; grid on
    % plot(pr_u1,pr_u2,'.r')
    % plot(sc_u1,sc_u2,'.c')
    % xlabel('u1'); ylabel('u2')
    
    
    % MSAC (robust variant of RANSAC)
    sampleSize  = 5;   % number of points to sample per trial
    maxDistance = 0.05; % max allowable distance for inliers
    
    fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
    % distance evaluation function
    evalLineFcn = @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
    
    points_pr = [pr_u1', pr_u2'];
    [model_pr, inlierIdx_pr] = ransac(points_pr,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
    yhat_pr = polyval(model_pr,pr_u1);
    
    points_sc = [sc_u1', sc_u2'];
    [model_sc, inlierIdx_sc] = ransac(points_sc,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
    yhat_sc = polyval(model_sc,sc_u1);
    
    % compute angle between fitted lines using normal to the lines
    n_pr = [1 model_pr(1)]; % normal vector = coefficients of x & y of line
    n_sc = [1 model_sc(1)]; % normal vector = coefficients of x & y of line
    costheta = dot(n_pr,n_sc)/(norm(n_pr)*norm(n_sc));
    theta = 180 - acosd(costheta); % incident angle between lines
    
    % subplot(212); hold on; grid on
    % plot(pr_u1,yhat_pr,'-r')
    % plot(sc_u1,yhat_sc,'-c')
    % legend('pr','sc','linefit on pr','linefit on sc')
    % 
    % Definitions
    slope = @(x,y) (y(end) - y(1))/(x(end) - x(1)); % slope of a line
    angle = @(m1,m2) atand((m1-m2)/(1+m1*m2));
    m1 = slope(pr_u1,yhat_pr); % slopoe of pre 
    m2 = slope(sc_u1,yhat_sc);
    theta2 = angle(m1,m2)
end

