% pnt_clust = pc_xyz(idxb,:);    % pnt_clust - cluster points
% pr_pnt = pnt_clust(1:idx,:);   % pr_pnt - preceeding points in the cluster
% sc_pnt = pnt_clust(idx:end,:); % sc_pnt - succeeding points in the cluster

num_pnts = size(pnt_clust,1);
idx = idx_mode;

pr_pnts = pnt_clust(1:idx,:);   % points preceeding the one having maxEigSum (including maxEigSum point)
sc_pnts = pnt_clust(idx:end,:);

[med,~ ] = mfames(pnt_clust'); % medoid
cent = mean(pc_xyz,1)'; % centroid
pnt_mode = pnt_clust(idx_mode,:);
const = pnt_mode';
pnt1 = pr_pnts - const';
pnt2 = sc_pnts - const';

% compute SVD
diff = pnt_clust-const';
Covar = diff'*diff;
Covar = Covar./num_pnts;
[U,E,~] = svd(Covar);
% [U,E,~] = svd((pnt_clust-const')');
u1 = U(:,1); u2 = U(:,2); u3 = U(:,3); % eigen vectors
e1 = E(1,1); e2 = E(2,2); e3 = E(3,3); % eigen values


[U_pr,sv_pr,~] = svd(pnt1');
[U_sc,sv_sc,~] = svd(pnt2');

u1_pr = U_pr(:,1);
u1_sc = U_sc(:,1);

% compute incident angle between the two lines
cosangle = @(a,b) dot(a,b)/(norm(a)*norm(b)); % cosine of incident angle

cosalpha = cosangle(u1_pr,u1_sc);
alpha = acosd(cosalpha)


% u-direction, v-direction
pr_u1 = u1'*pnt1'; % u1 component of pr
pr_u2 = u2'*pnt1'; % u2 component of pr
sc_u1 = u1'*pnt2'; % u1 component of sc
sc_u2 = u2'*pnt2'; % u2 component of sc


% MSAC (robust variant of RANSAC)
sampleSize  = 5;   % number of points to sample per trial
maxDistance = 0.5; % max allowable distance for inliers

fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
% distance evaluation function
evalLineFcn = @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);

rng(s);
% rng('default')

points_pr = [pr_u1', pr_u2'];
[model_pr, inlierIdx_pr] = ransac(points_pr,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
yhat_pr = polyval(model_pr,pr_u1);

points_sc = [sc_u1', sc_u2'];
[model_sc_old, inlierIdx_sc] = ransac(points_sc,fitLineFcn,evalLineFcn,sampleSize,maxDistance);
% [model_sc, ~ ] = linearfit(points_sc,inlierIdx_sc);
model_sc = polyfit(sc_u1(inlierIdx_sc)',sc_u2(inlierIdx_sc)',1);
% model_sc = robustfit(sc_u1(inlierIdx_sc)',sc_u2(inlierIdx_sc)','bisquare');
yhat_sc = polyval(model_sc,sc_u1);

% compute angle between fitted lines using normal to the lines
n_pr = [1 -model_pr(1)]; % normal vector = coefficients of x & y of line
n_sc = [1 -model_sc(1)]; % normal vector = coefficients of x & y of line
costheta = dot(n_pr,n_sc)/(norm(n_pr)*norm(n_sc));
theta = 180 - acosd(costheta) % incident angle between lines

vec_pr = [pr_u1(1), yhat_pr(1)] - [pr_u1(end), yhat_pr(end)];
vec_sc = [sc_u1(1), yhat_sc(1)] - [sc_u1(end), yhat_sc(end)];
vec_pr = vec_pr/norm(vec_pr);
vec_sc = vec_sc/norm(vec_sc);
dotprod = dot(vec_pr,vec_sc)


% n_sc2 = [1 -model_sc_robust(2)]; % normal vector = coefficients of x & y of line
% costheta2 = dot(n_pr,n_sc2)/(norm(n_pr)*norm(n_sc2));
% theta2 = 180 - acosd(costheta2) % incident angle between lines

% plot preceeding & succedding points
figure(501); clf; 
subplot(211); hold on; grid on
plot3(pnt1(:,1),pnt1(:,2),pnt1(:,3),'.r')
plot3(pnt2(:,1),pnt2(:,2),pnt2(:,3),'.c')
% axis equal

% plot eigen vectors
plot3([0 u1(1)],[0 u1(2)],[0 u1(3)],'-m')
plot3([0 u2(1)],[0 u2(2)],[0 u2(3)],'-g')
xlabel('x'); ylabel('y'); zlabel('z')
legend('pr','sc','u1','u2')
% axis equal

subplot(212); hold on; grid on
plot(pr_u1,pr_u2,'.r')
plot(sc_u1,sc_u2,'.c')

plot(pr_u1,yhat_pr,'-r')
plot(sc_u1,yhat_sc,'-c')
xlabel('u1'); ylabel('u2')
legend('pr','sc','pr linefit','sc linefit')

% %% Definitions
% slope = @(x,y) (y(end) - y(1))/(x(end) - x(1)); % slope of a line
% angle = @(m1,m2) atand((m1-m2)/(1+m1*m2));
% m1 = slope(pr_u1,yhat_pr); % slopoe of pre 
% m2 = slope(sc_u1,yhat_sc);
% theta2 = angle(m1,m2)

%%
% b1 = robustfit(pnt1(:,1:2),pnt1(:,3),'bisquare');
% b2 = robustfit(pnt2(:,1:2),pnt2(:,3),'bisquare');
% ct = cosangle(b1,b2);
% t = acosd(ct);

