function out_val = trim_scale(in_val,th)
%
% Trim by scaling
%

max_abs_val = max(abs(in_val(:)));
if max_abs_val > th
    out_val = in_val/max_abs_val*th;
else
    out_val = in_val;
end
