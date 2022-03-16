function chain = init_chain(varargin)
%
% Initialize a kinematic chain
%

% Parse options
ps = inputParser;
addParameter(ps,'name','kinematic_chain');
addParameter(ps,'dt',0.01);
parse(ps,varargin{:});
dt      = ps.Results.dt;
name    = ps.Results.name;

% Initialize a kinematic chain
chain                   = struct();
chain.name              = name;
chain.dt                = dt;
chain.joint             = struct(...
    'name','',...                   % joint name
    'p',cv([0,0,0]),...             % joint position
    'R',eye(3,3),...                % joint orientation
    'a',cv([0,0,0]),....            % rotation axis
    'type','fixed',...              % rotation type
    'p_offset',cv([0,0,0]),...      % position offset (w.r.t. parent joint)
    'R_offset',eye(3,3),...         % rotation offset (w.r.t. parent joint)
    'q',0,...                       % joint position
    'dq',0,...                      % dq/dt
    'ddq',0,...                     % d^2q/dt^2
    'q_diff',0,...                  % joint position different (q_diff/dt = dq)
    'q_prev',0,...                  % previous joint position
    'v',cv([0,0,0]),...             % linear velocity in {W}
    'vo',cv([0,0,0]),...            % linear part of spatial velocity
    'w',cv([0,0,0]),...             % angular velocity
    'dvo',cv([0,0,0]),...           % linear spatial acceleration
    'dw',cv([0,0,0]),...            % angular acceleration
    'u',0,...                       % joint torque
    'ext_f',cv([0,0,0]),...         % external force
    'parent',[],...                 % parent joint
    'childs',[],...                 % child joint
    'link_idx',[],...               % connected link index
    'limit',[-inf,+inf]...          % joint limit
    );
chain.joint_names       = {};
chain.n_joint           = 0;
chain.rev_joint_names   = {};
chain.n_rev_joint       = 0;

