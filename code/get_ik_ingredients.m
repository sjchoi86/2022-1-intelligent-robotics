function [J_use,ik_err] = get_ik_ingredients(chain,varargin)
%
% Get inverse kinematics ingredients
%

% Parse options
ps = inputParser;
addParameter(ps,'joint_names_to_ctrl',chain.rev_joint_names);
addParameter(ps,'joint_idxs_to_ctrl',idxs_cell(chain.joint_names,chain.rev_joint_names));
addParameter(ps,'joint_name_trgt','EE');
addParameter(ps,'T_trgt_goal',pr2t('',''));
addParameter(ps,'IK_P',1);
addParameter(ps,'IK_R',1);
addParameter(ps,'p_err_weight',1.0);
addParameter(ps,'w_err_weight',1.0);
addParameter(ps,'ik_err_th',10.0);
addParameter(ps,'DISREGARD_UNINFLUENTIAL_JOINT',0);
parse(ps,varargin{:});
joint_names_to_ctrl = ps.Results.joint_names_to_ctrl;
joint_idxs_to_ctrl  = ps.Results.joint_idxs_to_ctrl;
joint_name_trgt     = ps.Results.joint_name_trgt;
T_trgt_goal         = ps.Results.T_trgt_goal;
IK_P                = ps.Results.IK_P;
IK_R                = ps.Results.IK_R;
p_err_weight        = ps.Results.p_err_weight;
w_err_weight        = ps.Results.w_err_weight;
ik_err_th           = ps.Results.ik_err_th;
DISREGARD_UNINFLUENTIAL_JOINT = ps.Results.DISREGARD_UNINFLUENTIAL_JOINT;

% Compute the current target joint position from 'joint_name_trgt'
p_trgt_curr = chain.joint(idx_cell(chain.joint_names,joint_name_trgt)).p;
R_trgt_curr = chain.joint(idx_cell(chain.joint_names,joint_name_trgt)).R;
% Get the list of indices from root joint to the target joint
joint_idxs_route = get_idx_route(chain,joint_name_trgt);

% For positional IK, disregard revolute joints whose distance is small
if IK_P && DISREGARD_UNINFLUENTIAL_JOINT
    rmv_idxs = []; 
    for i_idx = 1:length(joint_idxs_route)
        joint_idx = joint_idxs_route(i_idx); 
        dist = norm(chain.joint(joint_idx).p-p_trgt_curr);
        if dist < eps
            rmv_idxs = cat(1,rmv_idxs,i_idx);
        end
    end
    joint_idxs_route(rmv_idxs) = []; % exclude unecessary joints 
end

% Intersect 'joint_idxs_route' with 'joint_idxs_to_control' to get actual using indices
joint_idxs_use = intersect(joint_idxs_route,joint_idxs_to_ctrl);
n_use = length(joint_idxs_use);
% Compute the Jacobian matrix (2.74)
n_ctrl = length(joint_idxs_to_ctrl);
J = zeros(6,n_ctrl);
for i_idx = 1:n_use % along the joint route
    joint_idx = joint_idxs_use(i_idx);
    parent = chain.joint(joint_idx).parent;      % parent joint index
    p_joint_ctrl = chain.joint(joint_idx).p;     % joint position
    R_offset = chain.joint(joint_idx).R_offset;  % current rotation offset
    % Rotation axis in {W}
    if isempty(parent)
        R_parent = eye(3,3);
    else
        R_parent = chain.joint(parent).R;
    end
    a = R_parent * R_offset * chain.joint(joint_idx).a;
    % 'idx_append': which column to append
    joint_name_append = chain.joint_names{joint_idx};
    idx_append = idx_cell(joint_names_to_ctrl,joint_name_append);
    J(:,idx_append) = [...
        cv(cross(a',p_trgt_curr-p_joint_ctrl));... % position part
        cv(a)... % orientation part (simply rotation axis in {W})
        ];
end % for i_idx = 1:n_use % along the joint route

% Compute the error
[p_trgt_goal,R_trgt_goal] = t2pr(T_trgt_goal);
J_use = [];
ik_err = [];
if IK_P
    J_use = [J_use; J(1:3,:)];
    p_err = p_trgt_goal-p_trgt_curr;
    if ~isfield(chain,'sz')
        chain.sz = get_chain_sz(chain); % get chain size if not exists
    end
    p_err = p_err / max(chain.sz.xyz_len); % normalize IK position error 
    ik_err = [ik_err; p_err_weight*p_err];
end
if IK_R
    J_use = [J_use; J(4:6,:)];
    w_err = R_trgt_curr * r2w(R_trgt_curr' * R_trgt_goal);
    ik_err = [ik_err; w_err_weight*w_err];
end
ik_err = trim_scale(ik_err,ik_err_th); % clamp IK error
