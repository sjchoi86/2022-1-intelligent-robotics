function chain = fk_chain(chain,idx_to)
%
% Forward kinematics
%
if isempty(idx_to)
    idx_to = get_topmost_idx(chain); % start from the topmost joint
end

% Update p and R of joint
idx_fr = chain.joint(idx_to).parent;
if ~isempty(idx_fr)
    joint_fr = chain.joint(idx_fr);
    joint_to = chain.joint(idx_to);
    
    % update p
    chain.joint(idx_to).p = joint_fr.R*joint_to.p_offset + joint_fr.p;
    
    % update R
    if isfield(joint_to,'a') % with rotational axis
        q = joint_to.q;
        a = joint_to.a;
        chain.joint(idx_to).R = ...
            joint_fr.R*joint_to.R_offset*rodrigues(a,q);
    else % without rotational axis
        chain.joint(idx_to).R = joint_fr.R*joint_to.R_offset;
    end
end

% Recursive
for child = chain.joint(idx_to).childs
    chain = fk_chain(chain,child);
end
