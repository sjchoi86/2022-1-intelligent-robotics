function chain = fv_chain(chain,idx_to)
%
% Forward velocity (p.93)
%
% We need 'q_diff' and 'dt' to be pre-computed.
%

if isempty(idx_to)
    idx_to = get_topmost_idx(chain); % start from the topmost joint
end

idx_fr = chain.joint(idx_to).parent;
if ~isempty(idx_fr)
    joint_fr = chain.joint(idx_fr);
    joint_to = chain.joint(idx_to);
    
    % Compute joint directional and angular velocities
    % (Kajita, 3.61)
    chain.joint(idx_to).v = joint_fr.v + cross(joint_fr.w,joint_fr.R*joint_to.p_offset);
    % (Katija, 3.62)
    q_dot_to = joint_to.q_diff/chain.dt;
    chain.joint(idx_to).w = joint_fr.w + (joint_to.R * joint_to.a * q_dot_to);
    
    % Compute the com velocity of the link attached to the joint
    link_index = chain.joint(idx_to).link_idx;
    if ~isempty(link_index)
        com_bar = chain.link(link_index).com_bar; % local com offset 
        v = chain.joint(idx_to).v;
        w = chain.joint(idx_to).w;
        % Compute directional and angular velocities of the link (Kajita, 3.64)
        chain.link(link_index).v = v + cross(w, joint_to.R*com_bar); 
        chain.link(link_index).w = w;
    end
end

% Recursive
for child = chain.joint(idx_to).childs
    chain = fv_chain(chain,child);
end
