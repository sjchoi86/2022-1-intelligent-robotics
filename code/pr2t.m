function T = pr2t(p,R)
%
% Get T from p and R
%

if isempty(p)
    p = cv([0,0,0]);
end

if isempty(R)
    R = eye(3,3);
end

T = [R,reshape(p,[3,1]);...
    0,0,0,1];
