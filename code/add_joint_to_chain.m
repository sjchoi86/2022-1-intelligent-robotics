function chain = add_joint_to_chain(chain,varargin)
%
% Add a joint to a kinematic chain
%

% Parse options
ps = inputParser;
addParameter(ps,'name','');                 % joint name
addParameter(ps,'p',cv([0,0,0]));           % joint position (further updated via FK)
addParameter(ps,'R',eye(3,3));              % joint rotation (further updated via FK)
addParameter(ps,'a',cv([0,0,0]));           % local axis of rotation
addParameter(ps,'p_offset',cv([0,0,0]));    % position offset 
addParameter(ps,'R_offset',eye(3,3));       % rotation offset 
addParameter(ps,'q',0);                     % joint position
addParameter(ps,'dq',0);                    % joint velocity 
addParameter(ps,'ddq',0);                   % joint acceleration 
addParameter(ps,'parent_name','');          % parent joint name 
addParameter(ps,'limit',[-inf,+inf]);       % joint position limit
addParameter(ps,'joint_names','');          % (auto-updated) joint names 
parse(ps,varargin{:});
name        = ps.Results.name;
p           = ps.Results.p;
R           = ps.Results.R;
a           = ps.Results.a;
p_offset    = ps.Results.p_offset;
R_offset    = ps.Results.R_offset;
q           = ps.Results.q;
dq          = ps.Results.dq;
ddq         = ps.Results.ddq;
parent_name = ps.Results.parent_name;
limit       = ps.Results.limit;
joint_names = ps.Results.joint_names;

% Increase jointn umber
chain.n_joint = chain.n_joint + 1;

% Name
if isempty(name)
    name = sprintf('joint_%02d',chain.n_joint);
end

% Joint axis type (automatically determined)
if isequal(a,cv([0,0,0]))
    type = 'fixed';
else
    type = 'revolute';
end

% Handle position offset
if isempty(p_offset)
    parent_idx = idx_cell(chain.joint_names,parent_name);
    if isempty(parent_idx)
        p_offset = cv([0,0,0]);
    else
        parent_pos = chain.joint(parent_idx).p;
        p_offset = cv(p) - cv(parent_pos);
    end
end

% Add joint
chain.joint(chain.n_joint) = struct(...
    'name',name,...
    'p',cv(p),'R',R,'a',cv(a),'type',type,...
    'p_offset',cv(p_offset),'R_offset',R_offset,...
    'q',q,'dq',dq,'ddq',ddq,'q_diff',0,'q_prev',0,...
    'v',cv([0,0,0]),'vo',cv([0,0,0]),'w',cv([0,0,0]),...
    'dvo',cv([0,0,0]),'dw',cv([0,0,0]),...
    'u',0,'ext_f',cv([0,0,0]),...
    'parent',[],'childs',[],...
    'link_idx',[],'limit',limit...
    );

% Add joint name
if isempty(joint_names)
    chain.joint_names{chain.n_joint} = name;
else
    % if 'joint_names' is given, then use it. 
    chain.joint_names = joint_names;
end

% Add parent
if ~isempty(parent_name)
    % Add parent
    parent_idx = idx_cell(chain.joint_names,parent_name);
    chain.joint(chain.n_joint).parent = [parent_idx];
    % Add the current joint as childs of the parent joint
    childs_of_parent = chain.joint(parent_idx).childs;
    childs_of_parent = [childs_of_parent, chain.n_joint];
    chain.joint(parent_idx).childs = childs_of_parent;
end

% Add revolute joint (if necesarry)
if isequal(type,'revolute')
    chain.n_rev_joint = chain.n_rev_joint + 1;
    chain.rev_joint_names{chain.n_rev_joint} = name;
end

% FK
chain = fk_chain(chain,'');
