function T = p2t(p)
%
% Get T from p and R
%

T = [eye(3,3),reshape(p,[3,1]);...
    0,0,0,1];
