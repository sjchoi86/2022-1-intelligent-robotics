function rpy_deg = r2rpy(R)
%
% Get roll, pitch, and yaw [in degree] from a rotation matrix
%  alpha: yaw
%  beta: pitch
%  gamma: roll 
%

r = atan2(R(3,2),R(3,3));
p = atan2(-R(3,1),sqrt(R(3,2)*R(3,2)+R(3,3)*R(3,3)));
y = atan2(R(2,1),R(1,1));

rpy_deg = [r,p,y];
