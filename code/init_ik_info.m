function ik_info = init_ik_info(chain,varargin)
%
% Initialize IK information
%

D2R = pi/180;

% Parse options
ps = inputParser;
addParameter(ps,'joint_names_to_ctrl',chain.rev_joint_names);
addParameter(ps,'ik_err_th',0.5);
addParameter(ps,'dq_th',45*D2R);
addParameter(ps,'step_size',0.5);
addParameter(ps,'lambda_rate',0.01);
addParameter(ps,'lambda_min',0.0001);
addParameter(ps,'lambda_max',0.1);
parse(ps,varargin{:});
joint_names_to_ctrl = ps.Results.joint_names_to_ctrl;
ik_err_th           = ps.Results.ik_err_th;
dq_th               = ps.Results.dq_th;
step_size           = ps.Results.step_size;
lambda_rate         = ps.Results.lambda_rate;
lambda_min          = ps.Results.lambda_min;
lambda_max          = ps.Results.lambda_max;

% Initialize IK information
ik_info.joint_names_to_ctrl = joint_names_to_ctrl;
ik_info.joint_idxs_to_ctrl  = idxs_cell(chain.joint_names,joint_names_to_ctrl);
ik_info.chain               = chain;
ik_info.n_trgt              = 0;     % number of IK targets
ik_info.trgt_joint_names    = {};    % IK target joint names
ik_info.trgt_types          = {};    % IK target types ('IK_P', 'IK_R')
ik_info.trgt_weights        = {};    % IK target weights
ik_info.trgt_coords         = {};    % IK target coordinates

% IK hyperparameters
ik_info.ik_err_th           = ik_err_th;
ik_info.dq_th               = dq_th;
ik_info.step_size           = step_size;
ik_info.lambda_rate         = lambda_rate;
ik_info.lambda_min          = lambda_min;
ik_info.lambda_max          = lambda_max;