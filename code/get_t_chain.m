function T = get_t_chain(chain,joint_name)
%
% Get joint coordinates
%

joint_idx = idx_cell(chain.joint_names,joint_name);
p = chain.joint(joint_idx).p;
R = chain.joint(joint_idx).R;
T = pr2t(p,R);
