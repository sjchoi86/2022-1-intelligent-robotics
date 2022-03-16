function [dq,joint_names_to_ctrl,J_use,ik_err,det_J] = one_step_ik(chain,ik_info,varargin)
%
% One-step IK
%

% Parse options
D2R = pi/180;
ps = inputParser;
addParameter(ps,'CONSIDER_JOINT_LIMIT',0);
addParameter(ps,'UNIT_DQ_HEURISTIC',1);
addParameter(ps,'unit_dq_rad',2*D2R);
addParameter(ps,'DISREGARD_UNINFLUENTIAL_JOINT',0);
parse(ps,varargin{:});
CONSIDER_JOINT_LIMIT = ps.Results.CONSIDER_JOINT_LIMIT;
UNIT_DQ_HEURISTIC    = ps.Results.UNIT_DQ_HEURISTIC;
unit_dq_rad          = ps.Results.unit_dq_rad;
DISREGARD_UNINFLUENTIAL_JOINT = ps.Results.DISREGARD_UNINFLUENTIAL_JOINT;


if ik_info.n_trgt == 0
    joint_names_to_ctrl = ik_info.joint_names_to_ctrl;
    dq = zeros(length(joint_names_to_ctrl),1);
    J_use = '';
    ik_err = '';
    det_J = '';
    return;
end

% Get dq for IK
[dq,J_use,ik_err,det_J] = get_dq_from_ik_info(chain,ik_info,...
    'DISREGARD_UNINFLUENTIAL_JOINT',DISREGARD_UNINFLUENTIAL_JOINT);

if CONSIDER_JOINT_LIMIT
    % Pre-update robot with dq and check joint limits
    q = get_q_chain(chain,ik_info.joint_names_to_ctrl);
    q_check = q + dq;
    joint_names_to_ctrl_valid = {}; joint_idxs_to_ctrl_valid = [];
    valid_cnt = 0; invalid_cnt = 0;
    for i_idx = 1:length(ik_info.joint_names_to_ctrl)
        q_check_i = q_check(i_idx);
        joint_name_to_check_i = ik_info.joint_names_to_ctrl{i_idx};
        joint_idx_to_check_i = idx_cell(chain.joint_names,joint_name_to_check_i);
        limit_i = chain.joint(joint_idx_to_check_i).limit;
        if (limit_i(1) < q_check_i) && (q_check_i < limit_i(2)) % check limit
            valid_cnt = valid_cnt + 1;
            joint_names_to_ctrl_valid{valid_cnt} = joint_name_to_check_i;
            joint_idxs_to_ctrl_valid = cat(1,joint_idxs_to_ctrl_valid,joint_idx_to_check_i);
        else
            % Not valid
            invalid_cnt = invalid_cnt + 1;
        end
    end
    if invalid_cnt > 0
        % If some joints exceeds its joint limit, recompute dq
        ik_info_valid = ik_info;
        ik_info_valid.joint_names_to_ctrl = joint_names_to_ctrl_valid;
        ik_info_valid.joint_idxs_to_ctrl = joint_idxs_to_ctrl_valid;
        [dq,J_use,ik_err,det_J] = get_dq_from_ik_info(chain,ik_info_valid);
        joint_names_to_ctrl = ik_info_valid.joint_names_to_ctrl;
        % Update robot with considering the joint limit
    else
        % Update robot as it is
        joint_names_to_ctrl = ik_info.joint_names_to_ctrl;
    end
else
    % Update robot without considering joint limits
    joint_names_to_ctrl = ik_info.joint_names_to_ctrl;
end

% Simple hueristics (default)
if UNIT_DQ_HEURISTIC
    dq = dq / norm(dq);
    dq = trim_scale(dq,unit_dq_rad);
end
