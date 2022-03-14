function pc_cyl = generate_cylinder(radius, height, cyl_origin, R, circle_count, circle_pc_count)
    % cyl_origin - center of the circle at the base of the cylider
    % R - rotation matrix to align the cylinder using specific rotation
    z_direction = [0 0 1];
    scales = linspace(0,height,circle_count);
    pc_circle = generate_circle(radius,circle_pc_count);    
    pc_cyl = [];
    for s = scales
        new_circle = cyl_origin + pc_circle + s.*z_direction;        
        pc_cyl = [pc_cyl; new_circle];
    end
    pc_cyl = pc_cyl*R';
end