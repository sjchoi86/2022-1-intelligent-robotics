function p = get_p_chain(chain,joint_name)
%
% Get joint coordinates
%

joint_idx = idx_cell(chain.joint_names,joint_name);
p = chain.joint(joint_idx).p;
