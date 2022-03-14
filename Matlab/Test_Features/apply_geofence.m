function fenced_pc = apply_geofence(pc,fence)
    % all fences are parallel normal to either x,y or z directions
    % there are lower & upper fences for each of x,y,z directions
    
    
    if nargin > 1
        fence_x_l = fence(1,1); % lower bound, normal to x-direction 
        fence_x_u = fence(1,2); % upper bound, normal to x-direction  
        fence_y_l = fence(2,1); % lower bound, normal to y-direction  
        fence_y_u = fence(2,2); % upper bound, normal to y-direction  
        fence_z_l = fence(3,1); % lower bound, normal to z-direction  
        fence_z_u = fence(3,2); % upper bound, normal to z-direction
    else
        fence_x_l = -inf; % lower bound, normal to x-direction 
        fence_x_u = inf;  % upper bound, normal to x-direction  
        fence_y_l = -inf; % lower bound, normal to y-direction  
        fence_y_u = inf;  % upper bound, normal to y-direction  
        fence_z_l = -inf; % lower bound, normal to z-direction  
        fence_z_u = 2;    % upper bound, normal to z-direction
    end
    fenced_pc = pc;

    % fences normal to x-direction
    if fence_x_l ~= -inf
        flag = fenced_pc(:,1) < fence_x_l;
        fidx = find(flag);
        fenced_pc(fidx,:) = [];
    end
    if fence_x_u ~= inf
        flag = fenced_pc(:,1) > fence_x_u;
        fidx = find(flag);
        fenced_pc(fidx,:) = [];
    end

    % fences normal to y-direction
    if fence_y_l ~= -inf
        flag = fenced_pc(:,2) < fence_y_l;
        fidx = find(flag);
        fenced_pc(fidx,:) = [];
    end
    if fence_y_u ~= inf
        flag = fenced_pc(:,2) > fence_y_u;
        fidx = find(flag);
        fenced_pc(fidx,:) = [];
    end

    % fences normal to z-direction
    if fence_z_l ~= -inf
        flag = fenced_pc(:,3) < fence_z_l;
        fidx = find(flag);
        fenced_pc(fidx,:) = [];
    end
    if fence_z_u ~= inf
        flag = fenced_pc(:,3) > fence_z_u;
        fidx = find(flag);
        fenced_pc(fidx,:) = [];
    end
end

