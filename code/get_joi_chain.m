function joi = get_joi_chain(chain)
%
% Get joints of intererst information (by searching joint names starting with 'joi_' 
% For example, 'joi_rh' => joi.type{1} = 'rh'
%
joi.n = 0; joi.idxs = []; joi.types = {};
for j_idx = 1:chain.n_joint
    joint_name = chain.joint_names{j_idx};
    if isequal(joint_name(1:min(4,length(joint_name))),'joi_')
        joi_type = joint_name(5:end);
        joi.idxs = [joi.idxs, j_idx];
        joi.n = joi.n + 1; joi.types{joi.n} = joi_type;
    end
end
