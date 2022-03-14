function [dx,dy,dz] = compute_derivative(pc_xyz,arg_plot)
    % computes dx,dy,dz
    
    dx = pc_xyz(2:end,1) - pc_xyz(1:end-1,1);
    dy = pc_xyz(2:end,2) - pc_xyz(1:end-1,2);
    dz = pc_xyz(2:end,3) - pc_xyz(1:end-1,3);
    N  = size(pc_xyz,1);
    index = [1:1:(N-1)]';
    mrkr_sz  = 3;
    mrkr_clr = "b";

    if lower(arg_plot) == 'drawfig'
        figure(8); 
        subplot(311); hold on; grid on
        scatter(index,dx,mrkr_sz,mrkr_clr);
        ylabel('dx');

        subplot(312); hold on; grid on
        scatter(index,dy,mrkr_sz,mrkr_clr);
        ylabel('dy');

        subplot(313); hold on; grid on
        scatter(index,dz,mrkr_sz,mrkr_clr);
        ylabel('dz');
        xlabel('index')
        sgtitle('derivative components');
    end
end

