function [cap,T,joint_name,joint_index] = get_chain_capsule(chain,link_index)
%
% Get chain capsule and coordinates from the link index
%

link = chain.link(link_index);
cap = link.capsule;
joint_index = link.joint_idx;
joint_name = chain.joint_names{joint_index};
T = pr2t(chain.joint(joint_index).p,chain.joint(joint_index).R);
