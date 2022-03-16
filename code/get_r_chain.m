function R = get_r_chain(chain,joint_name)
%
% Get joint coordinates
%

joint_idx = idx_cell(chain.joint_names,joint_name);
R = chain.joint(joint_idx).R;
