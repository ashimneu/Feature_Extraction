function pc_sphere = generate_sphere(radius,center,n,circle_pc_count)
    % Generates point cloud for sphere by slicing the sphere into latitude
    % and generating points of a circle for each latitude except poles
    % which are represented by single points.
    % n - [odd integer, >= 5] number of latitudes circles & points (points are for poles)
    % center    - position of center of sphere

    z_direction  = [0 0 1];
    pc_northpole = z_direction*radius;  % north pole
    pc_southpole = -z_direction*radius; % south pole
    equator   = generate_circle(radius,circle_pc_count);
    lat_count = n - 3;  % 3 is subtracted for 2 poles & 1 equator
    lat_count = lat_count/2; % latitude count for northern/southern hemisphere
    
    latitudes = linspace(0,pi/2,lat_count); % northern latitude: from 0 to 90 N
    radiuses  = radius.*cos(latitudes); % southern latitude: from 0 to 90 S
    levels    = linspace(0,radius,lat_count);
    levels(1) = [];
    radiuses(1)   = []; % remove first scale value to avoid redundant circle at equator circle
    levels(end)   = [];
    radiuses(end) = []; % remove last value to avoid redundant points at north pole
    
    % generate points for latitudes lines on northern hemisphere
    pc_sphere_N = [];
    for idx = 1:size(levels,2)
        level = levels(idx);
        lat_radius = radiuses(idx);
        pc_circle  = generate_circle(lat_radius,circle_pc_count); % points on a circle's circumference   
        translated_pnt_circle = center + pc_circle + level.*z_direction; % move to a latitude specified the scale value
        pc_sphere_N = [pc_sphere_N; translated_pnt_circle];
    end
    % generate points for latitudes lines on southern hemisphere
    pc_sphere_S = [];
    for idx = 1:size(levels,2)
        level = -levels(idx);
        lat_radius = radiuses(idx);
        pc_circle  = generate_circle(lat_radius,circle_pc_count); % points on a circle's circumference   
        translated_pnt_circle = center + pc_circle + level.*z_direction; % move to a latitude specified the scale value
        pc_sphere_S = [pc_sphere_S; translated_pnt_circle];
    end    
    pc_sphere = [pc_northpole; pc_sphere_N; equator; pc_sphere_S; pc_southpole];
end