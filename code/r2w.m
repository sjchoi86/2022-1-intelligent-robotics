function w = r2w(R)
%
% Get angular velocity [in radian] vector from a rotation matrix
% - This operation corresponds to w = wedge(ln(R)), a.k.a. a log map.
% - p.36 in 'Introduction to Humanoid Robotics' by Kajita et al.
% 
% The reverse can be computed by 'rodrigues'
% ex)
% w = r2w(R);
% R = rodrigues(uv(w),norm(w));
%
el = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
norm_el = norm(el);
if norm_el > eps
    w = atan2(norm_el, trace(R)-1)/norm_el * el;
elseif R(1,1)>0 && R(2,2)>0 && R(3,3)>0
    w = [0, 0, 0]';
else
    w = pi/2*[R(1,1)+1; R(2,2)+1; R(3,3)+1];
end

