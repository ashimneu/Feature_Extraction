function plot_components(p,fignum,pc_xyz,lzid)
    if numel(lzid) > 1
        fprintf('\n Plot of range, theta values is drawn for laser %2.0f only.\n',lzid(1));
    end
    
    figure(fignum); clf;
    pnt_num = size(pc_xyz,1);
    subplot(311); hold on; grid on; ylabel('x') 
    scatter(1:pnt_num,pc_xyz(:,1),p.mrk_size,'b')
    
    subplot(312); hold on; grid on; ylabel('y') 
    scatter(1:pnt_num,pc_xyz(:,2),p.mrk_size,'r')

    subplot(313); hold on; grid on; ylabel('z')
    scatter(1:pnt_num,pc_xyz(:,3),p.mrk_size,'g')
    xlabel('index')    
    sgtitle({['components, laser',num2str(lzid(1))]})
end

