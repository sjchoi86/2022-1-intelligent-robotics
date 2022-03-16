function chain_model = add_joi_to_coman(chain_model,varargin)
%
% Add joints of interest (JOI) with 'box_added' for IK and visualizaiton
%

D2R = pi/180;
% Parse options
p = inputParser;
addParameter(p,'parent_name_rh','RWrj2');
addParameter(p,'p_offset_rh',cv([0,0,0]));
addParameter(p,'rpy_offset_rh',[180,0,90]*D2R);
addParameter(p,'parent_name_lh','LWrj2');
addParameter(p,'p_offset_lh',cv([0,0,0]));
addParameter(p,'rpy_offset_lh',[180,0,90]*D2R);
addParameter(p,'parent_name_rw','RWrj2');
addParameter(p,'p_offset_rw',cv([0,0,0]));
addParameter(p,'rpy_offset_rw',[0,0,0]*D2R);
addParameter(p,'parent_name_lw','LWrj2');
addParameter(p,'p_offset_lw',cv([0,0,0]));
addParameter(p,'rpy_offset_lw',[0,0,0]*D2R);
addParameter(p,'parent_name_re','RElbj');
addParameter(p,'p_offset_re',cv([0,0,0]));
addParameter(p,'rpy_offset_re',[0,0,0]*D2R);
addParameter(p,'parent_name_le','LElbj');
addParameter(p,'p_offset_le',cv([0,0,0]));
addParameter(p,'rpy_offset_le',[0,0,0]*D2R);
addParameter(p,'parent_name_rs','RShLat');
addParameter(p,'p_offset_rs',cv([0,0,0]));
addParameter(p,'rpy_offset_rs',[0,0,0]*D2R);
addParameter(p,'parent_name_ls','LShLat');
addParameter(p,'p_offset_ls',cv([0,0,0]));
addParameter(p,'rpy_offset_ls',[0,0,0]*D2R);
addParameter(p,'parent_name_torso','base_joint');
addParameter(p,'p_offset_torso',cv([0,0,0]));
addParameter(p,'rpy_offset_torso',[0,0,0]*D2R);
addParameter(p,'parent_name_rf','r_leg_ft_joint');
addParameter(p,'p_offset_rf',cv([0,0,0]));
addParameter(p,'rpy_offset_rf',[90,0,90]*D2R);
addParameter(p,'parent_name_lf','l_leg_ft_joint');
addParameter(p,'p_offset_lf',cv([0,0,0]));
addParameter(p,'rpy_offset_lf',[90,0,90]*D2R);
addParameter(p,'parent_name_ra','RAnkSag');
addParameter(p,'p_offset_ra',cv([0,0,0]));
addParameter(p,'rpy_offset_ra',[0,0,0]*D2R);
addParameter(p,'parent_name_la','LAnkSag');
addParameter(p,'p_offset_la',cv([0,0,0]));
addParameter(p,'rpy_offset_la',[0,0,0]*D2R);
addParameter(p,'parent_name_rk','RKneeSag');
addParameter(p,'p_offset_rk',cv([0,0,0]));
addParameter(p,'rpy_offset_rk',[0,0,0]*D2R);
addParameter(p,'parent_name_lk','LKneeSag');
addParameter(p,'p_offset_lk',cv([0,0,0]));
addParameter(p,'rpy_offset_lk',[0,0,0]*D2R);
addParameter(p,'parent_name_rp','RHipLat');
addParameter(p,'p_offset_rp',cv([0,0,0]));
addParameter(p,'rpy_offset_rp',[0,0,0]*D2R);
addParameter(p,'parent_name_lp','LHipLat');
addParameter(p,'p_offset_lp',cv([0,0,0]));
addParameter(p,'rpy_offset_lp',[0,0,0]*D2R);
addParameter(p,'parent_name_head','');
addParameter(p,'p_offset_head',cv([0,0,0]));
addParameter(p,'rpy_offset_head',[0,0,0]*D2R);
addParameter(p,'hand_size',0.2);
addParameter(p,'wrist_size',0.2);
addParameter(p,'elbow_size',0.2);
addParameter(p,'shoulder_size',0.2);
addParameter(p,'torso_size',0.3);
addParameter(p,'foot_size',0.2);
addParameter(p,'ankle_size',0.2);
addParameter(p,'knee_size',0.2);
addParameter(p,'pelvis_size',0.2);
addParameter(p,'head_size',0.2);
%% Following should be common to all different robots
addParameter(p,'alpha',0.5);
addParameter(p,'color_rh','r');
addParameter(p,'color_lh','b');
addParameter(p,'color_rw','r');
addParameter(p,'color_lw','b');
addParameter(p,'color_re','r');
addParameter(p,'color_le','b');
addParameter(p,'color_rs','r');
addParameter(p,'color_ls','b');
addParameter(p,'color_torso','c');
addParameter(p,'color_rf','r');
addParameter(p,'color_lf','b');
addParameter(p,'color_ra','r');
addParameter(p,'color_la','b');
addParameter(p,'color_rk','r');
addParameter(p,'color_lk','b');
addParameter(p,'color_rp','r');
addParameter(p,'color_lp','b');
addParameter(p,'color_head','y');
parse(p,varargin{:});
parent_name_rh = p.Results.parent_name_rh;
p_offset_rh = p.Results.p_offset_rh;
rpy_offset_rh = p.Results.rpy_offset_rh;
parent_name_lh = p.Results.parent_name_lh;
p_offset_lh = p.Results.p_offset_lh;
rpy_offset_lh = p.Results.rpy_offset_lh;
parent_name_rw = p.Results.parent_name_rw;
p_offset_rw = p.Results.p_offset_rw;
rpy_offset_rw = p.Results.rpy_offset_rw;
parent_name_lw = p.Results.parent_name_lw;
p_offset_lw = p.Results.p_offset_lw;
rpy_offset_lw = p.Results.rpy_offset_lw;
parent_name_re = p.Results.parent_name_re;
p_offset_re = p.Results.p_offset_re;
rpy_offset_re = p.Results.rpy_offset_re;
parent_name_le = p.Results.parent_name_le;
p_offset_le = p.Results.p_offset_le;
rpy_offset_le = p.Results.rpy_offset_le;
parent_name_rs = p.Results.parent_name_rs;
p_offset_rs = p.Results.p_offset_rs;
rpy_offset_rs = p.Results.rpy_offset_rs;
parent_name_ls = p.Results.parent_name_ls;
p_offset_ls = p.Results.p_offset_ls;
rpy_offset_ls = p.Results.rpy_offset_ls;
parent_name_torso = p.Results.parent_name_torso;
p_offset_torso = p.Results.p_offset_torso;
rpy_offset_torso = p.Results.rpy_offset_torso;
parent_name_rf = p.Results.parent_name_rf;
p_offset_rf = p.Results.p_offset_rf;
rpy_offset_rf = p.Results.rpy_offset_rf;
parent_name_lf = p.Results.parent_name_lf;
p_offset_lf = p.Results.p_offset_lf;
rpy_offset_lf = p.Results.rpy_offset_lf;
parent_name_ra = p.Results.parent_name_ra;
p_offset_ra = p.Results.p_offset_ra;
rpy_offset_ra = p.Results.rpy_offset_ra;
parent_name_la = p.Results.parent_name_la;
p_offset_la = p.Results.p_offset_la;
rpy_offset_la = p.Results.rpy_offset_la;
parent_name_rk = p.Results.parent_name_rk;
p_offset_rk = p.Results.p_offset_rk;
rpy_offset_rk = p.Results.rpy_offset_rk;
parent_name_lk = p.Results.parent_name_lk;
p_offset_lk = p.Results.p_offset_lk;
rpy_offset_lk = p.Results.rpy_offset_lk;
parent_name_rp = p.Results.parent_name_rp;
p_offset_rp = p.Results.p_offset_rp;
rpy_offset_rp = p.Results.rpy_offset_rp;
parent_name_lp = p.Results.parent_name_lp;
p_offset_lp = p.Results.p_offset_lp;
rpy_offset_lp = p.Results.rpy_offset_lp;
parent_name_head = p.Results.parent_name_head;
p_offset_head = p.Results.p_offset_head;
rpy_offset_head = p.Results.rpy_offset_head;
alpha = p.Results.alpha;
color_rh = p.Results.color_rh;
color_lh = p.Results.color_lh;
color_rw = p.Results.color_rw;
color_lw = p.Results.color_lw;
color_re = p.Results.color_re;
color_le = p.Results.color_le;
color_rs = p.Results.color_rs;
color_ls = p.Results.color_ls;
color_torso = p.Results.color_torso;
color_rf = p.Results.color_rf;
color_lf = p.Results.color_lf;
color_ra = p.Results.color_ra;
color_la = p.Results.color_la;
color_rk = p.Results.color_rk;
color_lk = p.Results.color_lk;
color_rp = p.Results.color_rp;
color_lp = p.Results.color_lp;
color_head = p.Results.color_head;
hand_size = p.Results.hand_size;
wrist_size = p.Results.wrist_size;
elbow_size = p.Results.elbow_size;
shoulder_size = p.Results.shoulder_size;
torso_size = p.Results.torso_size;
foot_size = p.Results.foot_size;
ankle_size = p.Results.ankle_size;
knee_size = p.Results.knee_size;
pelvis_size = p.Results.pelvis_size;
head_size = p.Results.head_size;

% Right hand
jname2add_rh = 'joi_rh';
jname2add_rh_thumb = 'aux_right_thumb';
% Left hand
jname2add_lh = 'joi_lh';
jname2add_lh_thumb = 'aux_left_thumb';
% Right wrist
jname2add_rw = 'joi_rw';
% Left wrist
jname2add_lw = 'joi_lw';
% Right elbow
jname2add_re = 'joi_re';
% Left elbow
jname2add_le = 'joi_le';
% Right shoulder
jname2add_rs = 'joi_rs';
% Left shoulder
jname2add_ls = 'joi_ls';
% Torso
jname2add_torso = 'joi_torso';
jname2add_belly = 'aux_belly';
% Right foot
jname2add_rf = 'joi_rf';
jname2add_rf_ankle = 'aux_rf_ankle';
% Left foot
jname2add_lf = 'joi_lf';
jname2add_lf_ankle = 'aux_lf_ankle';
% Right ankle
jname2add_ra = 'joi_ra';
% Left ankle
jname2add_la = 'joi_la';
% Right knee
jname2add_rk = 'joi_rk';
% Left knee
jname2add_lk = 'joi_lk';
% Right pelvis
jname2add_rp = 'joi_rp';
% Left pelvis
jname2add_lp = 'joi_lp';
% Head
jname2add_head = 'joi_head';
jname2add_nose = 'aux_nose';

% Add joints and boxes 
% Add right hand
if ~isempty(parent_name_rh)
    % Add right palm
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rh,'parent_name',parent_name_rh,...
        'a',cv([0,0,0]),'p_offset',p_offset_rh,'R_offset',rpy2r(rpy_offset_rh));
    box_added = struct('xyz_min',hand_size*[-0.1,-0.3,0],'xyz_len',hand_size*[0.2,0.6,1.0],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_rh,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rh_palm_box',...
        'joint_name',jname2add_rh,'box_added',box_added);
    % Add right thumb
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rh_thumb,'parent_name',jname2add_rh,...
        'a',cv([0,0,0]),'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]));
    box_added = struct('xyz_min',hand_size*[-0.1,-0.25,0],'xyz_len',hand_size*[0.2,1.0,0.2],...
        'p_offset',hand_size*cv([0,0,0.15]),'R_offset',rpy2r([20,0,0]*D2R),...
        'color',color_rh,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rh_thumb_box',...
        'joint_name',jname2add_rh_thumb,'box_added',box_added);
end

% Add left hand (blue)
if ~isempty(parent_name_lh)
    % Add left palm
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_lh,'parent_name',parent_name_lh,...
        'a',cv([0,0,0]),'p_offset',p_offset_lh,'R_offset',rpy2r(rpy_offset_lh));
    box_added = struct('xyz_min',hand_size*[-0.1,-0.3,0],'xyz_len',hand_size*[0.2,0.6,1.0],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_lh,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rh_palm_box',...
        'joint_name',jname2add_lh,'box_added',box_added);
    % Add left thumb
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_lh_thumb,'parent_name',jname2add_lh,...
        'a',cv([0,0,0]),'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]));
    box_added = struct('xyz_min',hand_size*[-0.1,-0.25,0],'xyz_len',hand_size*[0.2,1.0,0.2],...
        'p_offset',hand_size*cv([0,0,0.15]),'R_offset',rpy2r([20,0,0]*D2R),...
        'color',color_lh,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rh_thumb_box',...
        'joint_name',jname2add_lh_thumb,'box_added',box_added);
end


% Add right wrist
if ~isempty(parent_name_rw)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rw,'parent_name',parent_name_rw,...
        'a',cv([0,0,0]),'p_offset',p_offset_rw,'R_offset',rpy2r(rpy_offset_rw));
    box_added = struct('xyz_min',wrist_size*[-0.1,-0.1,-0.1],'xyz_len',wrist_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_rw,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rw_box',...
        'joint_name',jname2add_rw,'box_added',box_added);
end

% Add left wrist
if ~isempty(parent_name_lw)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_lw,'parent_name',parent_name_lw,...
        'a',cv([0,0,0]),'p_offset',p_offset_lw,'R_offset',rpy2r(rpy_offset_lw));
    box_added = struct('xyz_min',wrist_size*[-0.1,-0.1,-0.1],'xyz_len',wrist_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_lw,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','lw_box',...
        'joint_name',jname2add_lw,'box_added',box_added);
end

% Add right foot (red)
if ~isempty(parent_name_rf)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rf,'parent_name',parent_name_rf,...
        'a',cv([0,0,0]),'p_offset',p_offset_rf,'R_offset',rpy2r(rpy_offset_rf));
    box_added = struct('xyz_min',foot_size*[-0.25,-0.1,-0.3],'xyz_len',foot_size*[0.5,0.2,1.0],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_rf,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rf_box',...
        'joint_name',jname2add_rf,'box_added',box_added);
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rf_ankle,'parent_name',jname2add_rf,...
        'a',cv([0,0,0]),'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]));
    box_added = struct('xyz_min',foot_size*[-0.15,-0.1,-0.15],'xyz_len',foot_size*[0.3,0.4,0.3],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_rf,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rf_ankle_box',...
        'joint_name',jname2add_rf_ankle,'box_added',box_added);
end

% Add left foot (blue)
if ~isempty(parent_name_lf)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_lf,'parent_name',parent_name_lf,...
        'a',cv([0,0,0]),'p_offset',p_offset_lf,'R_offset',rpy2r(rpy_offset_lf));
    box_added = struct('xyz_min',foot_size*[-0.25,-0.1,-0.3],'xyz_len',foot_size*[0.5,0.2,1.0],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_lf,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rf_box',...
        'joint_name',jname2add_lf,'box_added',box_added);
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_lf_ankle,'parent_name',jname2add_lf,...
        'a',cv([0,0,0]),'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]));
    box_added = struct('xyz_min',foot_size*[-0.15,-0.1,-0.15],'xyz_len',foot_size*[0.3,0.4,0.3],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_lf,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','rf_ankle_box',...
        'joint_name',jname2add_lf_ankle,'box_added',box_added);
end

% Add torso box
if ~isempty(parent_name_torso)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_torso,'parent_name',parent_name_torso,...
        'a',cv([0,0,0]),'p_offset',p_offset_torso,'R_offset',rpy2r(rpy_offset_torso));
    box_added = struct('xyz_min',torso_size*[-0.2,-0.3,-0.1],'xyz_len',torso_size*[0.4,0.6,0.9],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_torso,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','torso_box',...
        'joint_name',jname2add_torso,'box_added',box_added);
end
% Add torso belly
if ~isempty(parent_name_torso)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_belly,'parent_name',jname2add_torso,...
        'a',cv([0,0,0]),'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]));
    box_added = struct('xyz_min',torso_size*[-0.1,-0.4,-0.1],'xyz_len',torso_size*[0.4,0.8,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_torso,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,'name','torso_box2',...
        'joint_name',jname2add_belly,'box_added',box_added);
end

% Add right elbow
if ~isempty(parent_name_re)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_re,'parent_name',parent_name_re,...
        'a',cv([0,0,0]),'p_offset',p_offset_re,'R_offset',rpy2r(rpy_offset_re));
    box_added = struct('xyz_min',elbow_size*[-0.1,-0.1,-0.1],'xyz_len',elbow_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_re,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_re,'box_added',box_added,'name','re_box');
end
% Add left elbow
if ~isempty(parent_name_le)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_le,'parent_name',parent_name_le,...
        'a',cv([0,0,0]),'p_offset',p_offset_le,'R_offset',rpy2r(rpy_offset_le));
    box_added = struct('xyz_min',elbow_size*[-0.1,-0.1,-0.1],'xyz_len',elbow_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_le,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_le,'box_added',box_added,'name','le_box');
end

% Add right shoulder
if ~isempty(parent_name_rs)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rs,'parent_name',parent_name_rs,...
        'a',cv([0,0,0]),'p_offset',p_offset_rs,'R_offset',rpy2r(rpy_offset_rs));
    box_added = struct('xyz_min',shoulder_size*[-0.1,-0.1,-0.1],...
        'xyz_len',shoulder_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_rs,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_rs,'box_added',box_added,'name','rs_box');
end
% Add left shoulder
if ~isempty(parent_name_ls)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_ls,'parent_name',parent_name_ls,...
        'a',cv([0,0,0]),'p_offset',p_offset_ls,'R_offset',rpy2r(rpy_offset_ls));
    box_added = struct('xyz_min',shoulder_size*[-0.1,-0.1,-0.1],...
        'xyz_len',shoulder_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_ls,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_ls,'box_added',box_added,'name','ls_box');
end

% Add right pelvis
if ~isempty(parent_name_rp)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rp,'parent_name',parent_name_rp,...
        'a',cv([0,0,0]),'p_offset',p_offset_rp,'R_offset',rpy2r(rpy_offset_rp));
    box_added = struct('xyz_min',pelvis_size*[-0.1,-0.1,-0.1],...
        'xyz_len',pelvis_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_rp,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_rp,'box_added',box_added,'name','rp_box');
end
% Add left pelvis
if ~isempty(parent_name_lp)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_lp,'parent_name',parent_name_lp,...
        'a',cv([0,0,0]),'p_offset',p_offset_lp,'R_offset',rpy2r(rpy_offset_lp));
    box_added = struct('xyz_min',pelvis_size*[-0.1,-0.1,-0.1],...
        'xyz_len',pelvis_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_lp,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_lp,'box_added',box_added,'name','lp_box');
end

% Add right knee
if ~isempty(parent_name_rk)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_rk,'parent_name',parent_name_rk,...
        'a',cv([0,0,0]),'p_offset',p_offset_rk,'R_offset',rpy2r(rpy_offset_rk));
    box_added = struct('xyz_min',knee_size*[-0.1,-0.1,-0.1],'xyz_len',knee_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_rk,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_rk,'box_added',box_added,'name','rk_box');
end
% Add left knee
if ~isempty(parent_name_lk)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_lk,'parent_name',parent_name_lk,...
        'a',cv([0,0,0]),'p_offset',p_offset_lk,'R_offset',rpy2r(rpy_offset_lk));
    box_added = struct('xyz_min',knee_size*[-0.1,-0.1,-0.1],'xyz_len',knee_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_lk,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_lk,'box_added',box_added,'name','lk_box');
end

% Add right ankle
if ~isempty(parent_name_ra)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_ra,'parent_name',parent_name_ra,...
        'a',cv([0,0,0]),'p_offset',p_offset_ra,'R_offset',rpy2r(rpy_offset_ra));
    box_added = struct('xyz_min',ankle_size*[-0.1,-0.1,-0.1],'xyz_len',ankle_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_ra,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_ra,'box_added',box_added,'name','ra_box');
end
% Add left ankle
if ~isempty(parent_name_la)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_la,'parent_name',parent_name_la,...
        'a',cv([0,0,0]),'p_offset',p_offset_la,'R_offset',rpy2r(rpy_offset_la));
    box_added = struct('xyz_min',ankle_size*[-0.1,-0.1,-0.1],'xyz_len',ankle_size*[0.2,0.2,0.2],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_la,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_la,'box_added',box_added,'name','la_box');
end

% Add head
if ~isempty(parent_name_head)
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_head,'parent_name',parent_name_head,...
        'a',cv([0,0,0]),'p_offset',p_offset_head,'R_offset',rpy2r(rpy_offset_head));
    box_added = struct('xyz_min',head_size*[-0.35,-0.35,-0.1],'xyz_len',head_size*[0.7,0.7,0.8],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_head,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_head,'box_added',box_added,'name','head_box');
    chain_model = add_joint_to_chain(chain_model,...
        'name',jname2add_nose,'parent_name',parent_name_head,...
        'a',cv([0,0,0]),'p_offset',p_offset_head,'R_offset',rpy2r(rpy_offset_head));
    box_added = struct('xyz_min',head_size*[0.35,-0.1,0.25],'xyz_len',head_size*[0.2,0.2,0.15],...
        'p_offset',cv([0,0,0]),'R_offset',rpy2r([0,0,0]*D2R),...
        'color',color_head,'alpha',alpha,'ec','k');
    chain_model = add_link_to_chain(chain_model,...
        'joint_name',jname2add_nose,'box_added',box_added,'name','nose_box');
end


