clear; clc;

syms a b;

A = rot('x',a);
B = rot('z',b);

A*B

B*A

function rotation_matrix = rot(dim,alpha)
% alpha = deg2rad(alpha);
switch lower(dim)
    case 'x'
        rotation_matrix = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
    case 'y'
        rotation_matrix = [cos(alpha) 0 -sin(alpha); 0 1 0; sin(alpha) 0 cos(alpha)];
    case 'z'
        rotation_matrix = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
end
end
