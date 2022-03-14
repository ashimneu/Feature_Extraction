function [cor_pc, ln_pc] = scanline_feature_extract(PointCloud)

valid_range = 5; % threshold of valid slide window
sz_window = 11; % Size of the sliding window of each make sure to be odd integer
hf_sz = (sz_window-1)/2;

pc =  PointCloud.Location;

num_pnts = size(pc,1);
if(sz_window>num_pnts)
    disp('window too larger or equal to the number of points')
end

% initialize empty array to store computed feature values
pc_sw = zeros(sz_window,3); % pc of sliding window
feature_pc = [pc, zeros(num_pnts,1)]; % default type as 0;

% i = 695, plane before corner
% i = 716, corner
% i = 723 , orthogonal plane after corner

for i = 1:num_pnts
    first_pnt_idx = i - hf_sz;  % first point index
    last_pnt_idx = i + hf_sz;   % last point index


    % This if-statement below considers the case at the begining or end of
    % a scan. Basically it considers a scan as a circular scan. E.g., If 
    % the point chosen is the first point of a scan, then half of the slide 
    % window will be at the begining of the scan, the other half is the 
    % end of the scan.
    if first_pnt_idx <= 0 % at head
        head_idx = num_pnts - abs(first_pnt_idx);
        head_pc = pc(head_idx:num_pnts,:);
        pc_sw = [head_pc; pc(1:i+hf_sz,:)];
    elseif last_pnt_idx > num_pnts % at tail
        tail_idx = last_pnt_idx - num_pnts;
        tail_pc = pc(1:tail_idx,:);
        pc_sw = [ pc(i-hf_sz:num_pnts,:);tail_pc];
    else
        pc_sw =  pc(i-hf_sz:i+hf_sz,:);
    end

    %         figure(2)
    %         scatter3(pc_sw(:,1),pc_sw(:,2),pc_sw(:,3),2,'filled')
    %         axis([-15 25 -2 6 -4 16])
    %
    %         pause(0.1)


    dis_sep = norm(pc_sw(1,:)-pc_sw(end,:));
    if(dis_sep<valid_range) % only do analyse when all points are in a valid small range of slide window
        %         pnt_ref = mean(pc_sw); % calculate centroid of the points in the window
        pnt_ref = pc(i,:);
        diff = pc_sw - pnt_ref;  % calculate the difference between all points in the slide window with the chosen point
        H = diff'*diff;
        [~,S,~] = svd(H);

        if S(2,2)/S(3,3)>10000
            if S(1,1)/S(2,2)>10000
                feature_pc(i,4) = 1; % line
            elseif S(1,1)/S(2,2) < 10
                feature_pc(i,4) = 2; % corner
            end
        end        
        
    end
    %     [coeff, score, latent] = pca(pc_sw)
end

ln_pc = feature_pc( feature_pc(:,4)==1,:);
cor_pc = feature_pc( feature_pc(:,4)==2,:);
end