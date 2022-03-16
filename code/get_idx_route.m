function idx_route = get_idx_route(chain,joint_name_trgt,varargin)
%
% Get the list of indices from base to the target joint
%

% Parse options
p = inputParser;
addParameter(p,'EXCLUDE_TARGET',false); 
parse(p,varargin{:});
EXCLUDE_TARGET = p.Results.EXCLUDE_TARGET;

idx_target = idx_cell(chain.joint_names,joint_name_trgt);
if EXCLUDE_TARGET
    idx_route = []; % start from target, which will be reversed 
else
    idx_route = [idx_target]; % start from target, which will be reversed 
end

idx_temp = idx_target;
while 1
    parent_idx = chain.joint(idx_temp).parent;
    if isempty(parent_idx)
        break;
    end
    idx_route = cat(1,idx_route,parent_idx); % append parent
    idx_temp = parent_idx;
end
idx_route = flipud(idx_route); % flip to get base node to the end-effector

