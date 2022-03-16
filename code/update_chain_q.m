function chain = update_chain_q(chain,joint_names,q,varargin)
%
% Update the position of the kinematic chain
% 
% q: joint angles in radian
%

% Parse options
p = inputParser;
addParameter(p,'IGNORE_LIMIT',0);
addParameter(p,'FK',1);
addParameter(p,'FV',1);
parse(p,varargin{:});
IGNORE_LIMIT    = p.Results.IGNORE_LIMIT;
FK              = p.Results.FK;
FV              = p.Results.FV;

for i_idx = 1:length(joint_names)
    joint_idx = idx_cell(chain.joint_names,joint_names{i_idx});
    if isfield(chain.joint(joint_idx),'limit')
        limit = chain.joint(joint_idx).limit;
    else
        limit = [-inf,+inf];
    end
    if isempty(limit)
        limit = [-2*pi,+2*pi];
    end
    if IGNORE_LIMIT
        chain.joint(joint_idx).q = q(i_idx); 
    else
        chain.joint(joint_idx).q = min(max(q(i_idx),limit(1)),limit(2)); 
    end
    % Compute the joint position difference (if necesarry)
    if isfield(chain.joint(joint_idx),'q_prev')
        q_diff = chain.joint(joint_idx).q - chain.joint(joint_idx).q_prev;
        q_diff = mod(q_diff+pi,2*pi)-pi; % make sure the diff is between -pi ~ +pi
        chain.joint(joint_idx).q_diff = q_diff;
        chain.joint(joint_idx).q_prev = chain.joint(joint_idx).q;
    end
end

if FK
    chain = fk_chain(chain,'');
end

if FV
    chain = fv_chain(chain,'');
end

