function q = get_q_chain(chain,joint_names)
%
% Get the position (radian) of the kinematic chain
%

n_joint = length(joint_names);
q = zeros(n_joint,1);
for i_idx = 1:n_joint
    q(i_idx) = chain.joint(idx_cell(chain.joint_names,joint_names{i_idx})).q;
end
