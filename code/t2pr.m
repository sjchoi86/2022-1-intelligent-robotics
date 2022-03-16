function [p,R] = t2pr(T)
%
% Get p and R from T
%

if isempty(T)
    p = cv([0,0,0]);
    R = eye(3,3);
else
    p = T([1,2,3],4);
    R = T([1,2,3],[1,2,3]);
end

