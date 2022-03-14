function fout = plot_scanline(fignum,pc_xyz,cluster_ids,idx_s,idx_g,pc_s,pc_g,features,lzid)
    % fignum - figure number
    % pc_xyz - [numx3] points in scanline of 1 laser
    % cluster_ids - [numx1] cluster number of each point in pc_xyz
    % idx_s - index of corner points in pc_xyz
    % idx_g - index of non-corner points in pc_xyz
    % idx_s - corner points of pc_xyz
    % idx_g - non-corner points of pc_xyz
    % features - struct of features computed for points in pc_xyz
    % lzid - (optional) [1x1] laser number of the Lidar laser which was used to
    %         obtain pc_xyz points.


    % features
    angle   = features.angle;
    dotprod = features.dotprod;
    
    % indices
    pc_index = 1:size(pc_xyz,1); % indices of points in pc_xyz array

    % Parameters: figure
    sz_mrk = 75; % marker size    

    % Creat figure
    if fignum == 0
        figure(); hold on; grid on;
    else
        fighandle = figure(fignum); hold on; grid on;
    end
    
    % plot all points
    plt_all = scatter3(pc_xyz(:,1),pc_xyz(:,2),pc_xyz(:,3),sz_mrk,'.','MarkerEdgeColor',[0,1,0]);
    plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',pc_index);
    plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cluster',cluster_ids);
    plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dotprod',dotprod);
    plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('angle',angle);
    
    % plot corner points
    plt_c = scatter3(pc_s(:,1),pc_s(:,2),pc_s(:,3),sz_mrk,'.','MarkerEdgeColor',[1,0,0]);
    plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idx_s);
    plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cluster',cluster_ids(idx_s));
    plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dotprod',dotprod(idx_s));
    plt_c.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('angle',angle(idx_s));
    
    % plot non-corner points
    plt_nc = scatter3(pc_g(:,1),pc_g(:,2),pc_g(:,3),sz_mrk,'.','MarkerEdgeColor',[.3,.75,.93]);
    plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('index',idx_g);
    plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cluster',cluster_ids(idx_g));
    plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('dotprod',dotprod(idx_g));
    plt_nc.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('angle',angle(idx_g));
    
    xlabel('x'); ylabel('y'); zlabel('z')
    legend('outliers','cor.','non cor.')
    
    % load_fig_preferences();
    view([20 30])
    axis equal
    xlim([-100 100])
    ylim([-50  30])
    zlim([-30 5])

    
    if nargin == 9
        % append Lidar laser number of the pc to dataTip box.
        laserid    = lzid.*ones(size(pc_xyz,1),1);
        laserid_c  = lzid.*ones(numel(idx_s),1);
        laserid_nc = lzid.*ones(numel(idx_g),1);
        plt_all.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('laser',laserid);
        plt_c.DataTipTemplate.DataTipRows(end+1)   = dataTipTextRow('laser',laserid_c);
        plt_nc.DataTipTemplate.DataTipRows(end+1)  = dataTipTextRow('laser',laserid_nc);
    end

    if nargout == 1
        fout = fighandle;
    end
end

