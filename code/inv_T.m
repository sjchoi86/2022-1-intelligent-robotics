function T_out = inv_T(T)
%
% Get the inverse of T
%
p = t2p(T);
R = t2r(T);
T_out = [R' -R'*p; zeros(1,3), 1];
