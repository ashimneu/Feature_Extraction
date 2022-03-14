function plot_range_theta(p,fignum,range,theta,lzid)
    % fig num - figure number
    if numel(lzid) > 1
        fprintf('\n Plot of range, theta values is drawn for laser %2.0f only.\n',lzid(1));
    end

    pnt_num = numel(range); % number of points in a scanline

    % Plot range & range        
    diff_range = range(2:end) - range(1:end-1);
    diff_theta = theta(2:end) - theta(1:end-1);
    gr = 2; gc = 2;
    
    
    figure(fignum); clf
    
    % Plot range 
    subplot(gr,gc,1)
    hold on; grid on
    scatter(1:pnt_num,range,p.mrk_size,'b')
    ylabel('range')
    
    % Plot range difference
    subplot(gr,gc,2)
    hold on; grid on
    scatter(1:pnt_num-1,diff_range,p.mrk_size,'b')
    ylabel('range diff.')

    % Plot theta 
    subplot(gr,gc,3)
    hold on; grid on
    scatter(1:pnt_num,theta,p.mrk_size,'b')
    xlabel('point index')
    ylabel('theta')
    
    % Plot theta difference
    subplot(gr,gc,4)
    hold on; grid on
    scatter(1:pnt_num-1,diff_theta,p.mrk_size,'b')
    xlabel('diff. index')
    ylabel('theta diff.')
    
    sgtitle({['scanline, laser',num2str(lzid(1))]})
end

