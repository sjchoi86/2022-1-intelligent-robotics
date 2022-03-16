function ik_info = add_ik_info(ik_info,varargin)
%
% Add IK information
%

% Parse options
ps = inputParser;
addParameter(ps,'joint_name','');
addParameter(ps,'type','IK_P');
addParameter(ps,'weight',1);
addParameter(ps,'coord',pr2t('',''));
parse(ps,varargin{:});
joint_name  = ps.Results.joint_name;
type        = ps.Results.type;
weight      = ps.Results.weight;
coord       = ps.Results.coord;

% Add information
ik_info.n_trgt = ik_info.n_trgt + 1;
ik_info.trgt_joint_names{ik_info.n_trgt} = joint_name;
ik_info.trgt_types{ik_info.n_trgt}       = type;
ik_info.trgt_weights{ik_info.n_trgt}     = weight;
ik_info.trgt_coords{ik_info.n_trgt}      = coord;
