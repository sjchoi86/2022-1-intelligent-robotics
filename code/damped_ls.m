function [dq,det_J] = damped_ls(J_use,ik_err,varargin)
%
% Damped least square
%

% Parse options
p = inputParser;
addParameter(p,'lambda_rate',0.01);
addParameter(p,'lambda_min',1e-6);
addParameter(p,'lambda_max',1e-0);
addParameter(p,'step_size',0.1);
addParameter(p,'dq_th',20*pi/180);
parse(p,varargin{:});
lambda_rate = p.Results.lambda_rate;
lambda_min  = p.Results.lambda_min;
lambda_max  = p.Results.lambda_max;
step_size   = p.Results.step_size;
dq_th       = p.Results.dq_th;

% Damped least square
ik_err_avg = mean(abs(ik_err));
lambda = lambda_rate*ik_err_avg + lambda_min; % damping term
n_ctrl = size(J_use,2);

% Lambda scheduling 
idx_nz = abs(sum(J_use,1))>0.1;
J_use_nz = J_use(:,idx_nz);
det_J = det(J_use_nz*J_use_nz');
if det_J > 1e-3
    lambda = 1e-6;
elseif det_J < 1e-20
    lambda = lambda_max;
end

dq_raw = (J_use'*J_use + lambda*eye(n_ctrl,n_ctrl)) \ J_use' * ik_err;
dq = trim_scale(step_size*dq_raw,dq_th);
