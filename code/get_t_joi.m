function T_joi = get_t_joi(chain,joi)
%
% Get T of chain's joints of interest
%

T_joi = struct();
for i_idx = 1:length(joi.types)
    joi_type    = joi.types{i_idx};
    joint_idx   = joi.idxs(i_idx);
    p_i = chain.joint(joint_idx).p;
    R_i = chain.joint(joint_idx).R;
    T_i = pr2t(p_i,R_i);
    T_joi  = setfield(T_joi,joi_type,T_i);
end
