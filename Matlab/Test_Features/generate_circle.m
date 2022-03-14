function pc_circle = generate_circle(radius,n)
    % n - number of points on the circle
    sector_angles = linspace(0,2*pi,n)';
    pc_circle = radius.*[cos(sector_angles), sin(sector_angles)];
    pc_circle = [pc_circle zeros(n,1)];
end