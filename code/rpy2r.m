function R = rpy2r(rpy_rad)
%
% Euler angles (roll, pitch, and yaw) in radian to a rotation matrix
%
r_rad = rpy_rad(1);
p_rad = rpy_rad(2);
y_rad = rpy_rad(3);

cos_r = cos(r_rad); sin_r = sin(r_rad);
cos_p = cos(p_rad); sin_p = sin(p_rad);
cos_y = cos(y_rad); sin_y = sin(y_rad);

R = [...
    cos_y*cos_p,   -sin_y*cos_r+cos_y*sin_p*sin_r, sin_y*sin_r+cos_y*sin_p*cos_r ; ...
    sin_y*cos_p,    cos_y*cos_r+sin_y*sin_p*sin_r,  -cos_y*sin_r+sin_y*sin_p*cos_r ; ...
    -sin_p,         cos_p*sin_r,                    cos_p*cos_r...
    ];
