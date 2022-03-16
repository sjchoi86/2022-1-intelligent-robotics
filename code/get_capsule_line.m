function cl = get_capsule_line(T_cap,cap)
%
% returns a line segment (two points) from a capsule handle
%
% calculate positions
% cl.p1 = reshape(p,[3,1]) - 0.5*height*R*u;
% cl.p2 = reshape(p,[3,1]) + 0.5*height*R*u;
%

T = T_cap*cap.T_offset;
p = t2p(T);
R = t2r(T);
height = cap.height;
u = [0 0 1]';

cl.p1 = reshape(p,[3,1]) + 0.5*height*R*u;
cl.p2 = reshape(p,[3,1]) - 0.5*height*R*u;
