clc; clear;
p.add_noise = false;
p.eb_fig_moments = false;

corner_angles = 45; %30:30:360; % 2dcorner angle
gaps = 0:1:10; % gap between two adjacent lines
test_variables = corner_angles;
for j = 1:numel(test_variables)
    test_var = test_variables(j);
pt_num = 100; % number of points on a line or a circle.
sz_window = 25;
pnt_separation = .4; % point separation
ratio = pnt_separation/sz_window; % ratio of distance between successive point in a scanline to sliding window length
fprintf('\nratio : %3.3f \n',ratio);


% 2D Geometical Shapes
% generate points on the circumference of a circle
pc_circle = generate_circle(pt_num,25);
% generate a line of points along x-axis. Columns denote x,y,z coordinates
pc_line   = generate_line(p,[0.5 0 0],[1 0 0],pt_num,pnt_separation);
% generate two line of points along x-axis with a gap (10 units) in between them.
line1 = generate_line(p,[0 0 0],[1 0 0],pt_num,1);
line2 = generate_line(p,line1(end,:)+[10 0 0],[1 0 0],pt_num,1); % disable this line if comparing results for various gap values
% line2 = generate_line(line1(end,:)+[test_var 0 0],[1 0 0],pt_num,1); % enable this line if comparing results for various gap values

pc_linewithgap  = [line1; line2];
% generate a 2d corner (two connecting lines)
% pc_corner2d = generate_2dcorner(pc_line,90); % disable this line if comparing results for various corner angles
pc_corner2d = generate_2dcorner(p,pc_line,test_var); % enable this line if comparing results for various corner angles

% 3D Geometrical Shapes
% generate a plane by repeating a line of points
pc_plane  = generate_plane(pc_line,10); 
% generate points resembling an edge by rotating pc_plane about Y-axis 
pc_edge1  = generate_edge(pc_plane,[0,90,0]); 
% generate points resembling an edge by translating a 2d corner pointcloud along a direction
pc_edge2  = generate_edge_from_2dcorner(pc_corner2d,[0,0,0],10,5,[0,0,0]);
% generate points that resemble corner by rotating pc_plane about Y-axis & X-axis 
pc_corner3d = generate_3dcorner(pc_plane,[0,90,0],[-90,0,0]); 
% generate cylinder points of given radius, height, origin, orientation rot
% matrix, count of circles & points on each circle.
R =  rot('y',90); % 
pc_cylinder = generate_cylinder(5, 50, [0 0 0], R, 10, 30);
% generate sphere points of given radius
pc_sphere = generate_sphere(8,[0 0 0],35,20);

% Choose PointCloud for feature computation
% pc = pc_linewithgap;
pc = pc_corner2d*rot('y',150)'; % 2nd moments based corner detections works for 90deg corners
% pc = pc_edge2; % 2nd moments based corner detections works for 90deg corners
% pc = pc_corner3d;

pr_pnts = pc(1:248,:);
sc_pnts = pc(249:end,:);

% PCA over preceeding points (step 2a)
[V_pr,E_pr,~] = compute_pca(pr_pnts');
E_pr1 = E_pr(1,1); E_pr2 = E_pr(2,2);
V_pr1 = V_pr(:,1); V_pr2 = V_pr(:,2);

% PCA over succeeding points (step 2b)
[V_sc,E_sc,~] = compute_pca(sc_pnts');
E_sc1 = E_sc(1,1); E_sc2 = E_sc(2,2);
V_sc1 = V_sc(:,1); V_sc2  = V_sc(:,2);

% (step 3) check if preceeding & succeeding points are linear points
E_threshold = 2; % singular value threshold to satisfy single dimensionality of the points
pr_ratio2 = (E_pr1 - E_pr2 )/E_pr2;
sc_ratio2 = (E_sc1 - E_sc2 )/E_sc2;
pr_is_linear2 = pr_ratio2 > E_threshold; % threshold check
sc_is_linear2 = sc_ratio2 > E_threshold; % threshold check
dotprod2 = abs(dot(V_pr1,V_sc1)); % dot product
cos_alpha = dot(V_pr1,V_sc1)/(norm(V_pr1)*norm(V_sc1));
alpha = acosd(cos_alpha); % angle between lines

fprintf('\nPCA: \n')
fprintf('E1_pr = %1.3f, E2_pr = %1.3f, rat2_pr = %3.2f \n',[E_pr1,E_pr2,pr_ratio2])
fprintf('E1_sc = %1.3f, E2_sc = %1.3f, rat2_sc = %3.2f \n',[E_sc1,E_sc2,sc_ratio2])
fprintf('dp = %1.3f \n',dotprod2)
fprintf('angle = %1.3f \n',alpha)
[pr_is_linear2,sc_is_linear2]
direction2 = [V_pr1';V_sc1']
direction2 = [V_pr1';V_sc1'];


% Apply Geofence
fence = [-inf(3,1) inf(3,1)];
% fence(2,2) = 0; 
pc = apply_geofence(pc,fence);
N  = size(pc,1);
centroid = sum(pc,1)/N;
diff  = pc - centroid;
Covar = diff'*diff./(N+1);

% compute eigen values & eigen vectors
[Evec, Eval, ~] = eig(Covar);
[Eval,Idx] = sort(diag(Eval),'descend'); 
eig1 = Eval(1); evec1 = Evec(:,Idx(1));
eig2 = Eval(2); evec2 = Evec(:,Idx(2));
eig3 = Eval(3); evec3 = Evec(:,Idx(3));

fprintf('\nEigen val. 1 : %3.3f \n',eig1);
fprintf('Eigen val. 2 : %3.3f \n',eig2);
fprintf('Eigen val. 3 : %3.3f \n',eig3);

% Compute Geometric Features (for points in 3D pc)
L = (eig1 - eig2)/eig1;
P = (eig2 - eig3)/eig1;
S = eig3/eig1;
O_3d    = (eig1*eig2*eig3)^(1/3);
A_3d    = (eig1-eig3)/eig1;
Sum_3d  = eig1+eig2+eig3;
SV_3d   = eig3/(Sum_3d); % surface variation in 3D
etrp_3d = -(eig1*log(eig1) + eig2*log(eig2) + eig3*log(eig3)); % Eigen Entropy

% Compute Adapted Geometric Features (for points in 2D pc)
C       = eig2/eig1;
O_2d    = (eig1*eig2)^(1/2);
A_2d    = (eig1-eig2)/eig1;
Sum_2d  = eig1+eig2;
SV_2d   = eig2/(Sum_2d); % surface variation in 2D

% Compute Moments
% FOFA - 1st order, 1st axis
% FOSA - 1st order, 2nd axis
% SOFA - 2nd order, 1st axis
% SOSA - 2nd order, 2nd axis
[FOFA,FOSA,SOFA,SOSA] = compute_moment_using_sw(pc,sz_window);

% Compute Smoothness
smoothness = compute_smoothness_using_sw(pc,sz_window);

% Compute Linearity
linearity = compute_linearity_using_sw(pc,sz_window);

% % Print Geometric Features (3D)
% fprintf('\n');
% fprintf('Geometric Features in 3D \n');
% fprintf('Linearity   : %3.3f \n',L);
% fprintf('Planarity   : %3.3f \n',P);
% fprintf('Sphericity  : %3.3f \n',S);
% fprintf('Omnivar.    : %3.3f \n',O_3d);
% fprintf('Anisotropy  : %3.3f \n',A_3d);
% fprintf('Eigval Sum  : %3.3f \n',Sum_3d);
% fprintf('Surf. Var.  : %3.3f \n',SV_3d);
% fprintf('EigEntropy  : %3.3f \n',etrp_3d);
% 
% % Print Geometric Features (2D)
% fprintf('\n');
% fprintf('Geometric Features in 2D \n');
% fprintf('Linearity   : %3.3f \n',L);
% fprintf('Circularity : %3.3f \n',C);
% fprintf('Omnivar.    : %3.3f \n',O_2d);
% fprintf('Anisotropy  : %3.3f \n',A_2d);
% fprintf('Eigval Sum  : %3.3f \n',Sum_2d);
% fprintf('Surf. Var.  : %3.3f \n',SV_2d);

% Plot the point cloud & eigen vectors
figure(1); clf; hold on; grid on
plot3(pc(:,1),pc(:,2),pc(:,3),'bo');
xlabel('x')
ylabel('y')
zlabel('z')


% plot x & y axis lines on Figure(1)
x_axis = linspace(-10,10,10);
y_axis = linspace(0,0,10);
plot(x_axis,y_axis,'k-|','HandleVisibility','off');
plot(y_axis,x_axis,'k-|','HandleVisibility','off');

% plot all eigen vectors
plot3([0,evec1(1)],[0,evec1(2)],[0,evec1(3)],'r-')
plot3([0,evec2(1)],[0,evec2(2)],[0,evec2(3)],'g-')
plot3([0,evec3(1)],[0,evec3(2)],[0,evec3(3)],'m-')
legend('pc','e\_vec 1','e\_vec 2','e\_vec 3')
view([20 30])
axis equal

if p.eb_fig_moments
    % Plot 1st & 2nd moments
    figure(2); clf; 
    subplot(221);hold on; grid on
    plot(1:N,FOFA,'.')
    ylabel('Moment')
    title('1st Moment, 1st axis')
    
    subplot(222);hold on; grid on
    plot(1:N,FOSA,'.')
    title('1st Moment, 2nd axis')
    
    subplot(223);hold on; grid on
    plot(1:N,SOFA,'.')
    xlabel('index');
    ylabel('Moment')
    title('2nd Moment, 1st axis')
    
    subplot(224);hold on; grid on
    plot(1:N,SOSA,'.')
    xlabel('index');
    title('2nd Moment, 2nd axis')
    axis equal
    
    % Plot Smoothness
    figure(3); clf; 
    subplot(211); hold on; grid on
    plot(1:N,smoothness,'.')
    ylabel('Moment')
    title('Smoothness')
    axis equal
    
    subplot(212); hold on; grid on
    plot(1:N,linearity,'.')
    ylabel('Linearity')
    xlabel('index');
    title('Linearity')
end
end

function pc_line = generate_line(p,origin,direction,length, pnt_separation)
    ADD_NOISE = p.add_noise;
    scale = 0:pnt_separation:length-1;
    num_scale = numel(scale);
    direction = direction/norm(direction); % normalize direction vector
    temp_direc = repmat(direction,num_scale,1);
    pc_line = origin + scale'.*temp_direc;
    
    if ADD_NOISE
        % generate noise
        num_pnt = size(pc_line,1);
        noise_std = diag([0.15, 0.15, 0.45]); % noise standard deviations
        noise = randn(num_pnt,3)*noise_std';   
        
        pc_line = pc_line + noise; % add noise
    end
end

function pc_plane = generate_plane(pc_line,n)
    % pc_line - line points parallel to x-axis
    % n - number of times pc_line is repeated

    N = size(pc_line,1);
    
    pc_plane = [pc_line];
    for i = 1:n-1
        new_line = pc_line + [zeros(N,1), i.*ones(N,1), zeros(N,1)];
        pc_plane = [pc_plane; new_line];
    end    
end

function pc_output = generate_2dcorner(p,pc_2dline,gamma)
    % gamma - [degree] rotation angles about z direction
    % pc_2dline - [Nx3] pc of line in x-y plane where z entries of all
    % points are 0.
    
    ADD_NOISE = p.add_noise;

    N = size(pc_2dline,1);
    
    % make z-entries of all points 0 even if they are not.
    pc_2dline(:,3) = zeros(N,1);
    
    % rotate the line about desired angle about z-axis
    R = rot('z',gamma);    
    pc_line_rotated = pc_2dline*R';    

    % flip entries for sequential indexing. Sequential indexing of points is
    % required for computation of moments.
    pc_line_rotated = flip(pc_line_rotated,1);
    
    pc_corner_2d = [pc_line_rotated; pc_2dline];
    pc_output    = pc_corner_2d;
    
    if ADD_NOISE
        % generate noise
        num_pnt = size(pc_corner_2d,1);
        noise_std = diag([0.15, 0.15, 0.45]); % noise standard deviations
        noise = randn(num_pnt,3)*noise_std';   
        % add noise
        pc_output = pc_corner_2d + noise; 
    end
end

function pc_edge = generate_edge(pc_plane,abg)
    % generate edge points by rotating a plane pc about x,y,z axis 
    % abg - alpha, beta, gamma - rotation angles about x,y,z
    alpha = abg(1); 
    beta  = abg(2); 
    gamma = abg(3);
    R = rot('x',alpha)*rot('y',beta)*rot('z',gamma);    
    pc_plane_rotated = pc_plane*R';
    
    % flip entries for sequential indexing. Sequential indexing of points is
    % required for computation of moments.
    pc_plane_rotated = flip(pc_plane_rotated,1);
    
    pc_edge = [pc_plane_rotated; pc_plane];
end

function pc_2dcorner = generate_edge_from_2dcorner(pc_corner_2d,origin,length,rep_num,abg)
    % generate edge point cloud using translated copy of given pc of 2d
    % corner along a direction vector.
    
    % pc_corner_2d - (Nx3) pc of corner points lying in a plane.
    % Points in the array need to be sequentially ordered corresponding 
    % to their spatial position during the scan.
    % length - desired length of the edge
    % rep_num - number of repetition of pc_corner_2d
    % abg - alpha, beta, gamma - rotation angles about x,y,z
    % abg should be provided such that the resultant rotation matrix 
    % aligns pc_corner_2d with the desired direction of the edge.
      
    scales = linspace(0,length,rep_num);
    direction = [0 0 1];
    pc_temp = [];
    for i = 1:rep_num
        scale = scales(i);
        pc_translated =  origin + scale.*direction + pc_corner_2d;
        pc_temp = [pc_temp; pc_translated];
    end
    alpha = abg(1); 
    beta  = abg(2); 
    gamma = abg(3);
    R = rot('x',alpha)*rot('y',beta)*rot('z',gamma); 
    pc_rotated = pc_temp*R';
    
    pc_2dcorner = pc_rotated ; % + noise; % noise is added to the generated clean pc
end

function pc_corner = generate_3dcorner(pc_plane,abg1,abg2)
    % abg - alpha, beta, gamma - rotation angles about x,y,z
    alpha = abg1(1); 
    beta  = abg1(2); 
    gamma = abg1(3);
    R1 = rot('x',alpha)*rot('y',beta)*rot('z',gamma);
    pc_plane_rotated_1 = R1*pc_plane';

    alpha = abg2(1); 
    beta  = abg2(2); 
    gamma = abg2(3);
    R2 = rot('x',alpha)*rot('y',beta)*rot('z',gamma);
    pc_plane_rotated_2 = R2*pc_plane';

    pc_corner = [pc_plane; pc_plane_rotated_1'; pc_plane_rotated_2'];
end


