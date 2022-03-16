function out = uv(vec)
%
% To a column

nv = norm(vec);
if nv < 1e-8
    out = 0*vec;
else
    out = vec / nv;
end
