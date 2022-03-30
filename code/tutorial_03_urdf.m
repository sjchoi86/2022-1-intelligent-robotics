ccc
%% How to parse a URDF file
ccc
% Parse raw xml file
robot_name = 'iiwa7';
urdf_path = sprintf('../urdf/%s/%s_urdf.xml',robot_name,robot_name);
s = xml2str(urdf_path);
robot = s.robot;

% Parse basic joint information
chain = init_chain('name',robot_name);
n_joint = length(robot.joint); % number of joints

% Basic joint and link information
joint_names       = cell(1,n_joint);
parent_link_names = cell(1,n_joint);
child_link_names  = cell(1,n_joint);
rev_joint_names = {}; n_rev_joint = 0;
for i_idx = 1:n_joint % for all joints
    joint_i = robot.joint{i_idx};
    joint_name = joint_i.Attributes.name;
    joint_type = joint_i.Attributes.type;
    joint_names{i_idx} = joint_name;
    parent_link_names{i_idx} = joint_i.parent.Attributes.link;
    child_link_names{i_idx} = joint_i.child.Attributes.link;
    if isequal(joint_type,'revolute')
        n_rev_joint = n_rev_joint + 1;
        rev_joint_names{n_rev_joint} = joint_name;
    end
end % for i_idx = 1:n_joint % for all joints
chain.joint_names     = joint_names;
chain.n_joint         = n_joint;
chain.rev_joint_names = rev_joint_names;
chain.n_rev_joint     = n_rev_joint;

% Parse detailed joint information
for i_idx = 1:n_joint % for all joints
    joint_i = robot.joint{i_idx};
    joint_name = joint_i.Attributes.name;
    joint_type = joint_i.Attributes.type;
    % Joint position limit
    if isfield(joint_i,'limit')
        limit_lower = joint_i.limit.Attributes.lower;
        upper = joint_i.limit.Attributes.upper;
        if limit_lower(1) == '$', limit_lower = eval(limit_lower(3:end-1));
        else, limit_lower = eval(limit_lower);
        end
        if upper(1) == '$', upper = eval(upper(3:end-1));
        else, upper = eval(upper);
        end
        limit = [limit_lower,upper];
    else
        if isequal(joint_type,'revolute')
            limit = [-inf,+inf]; % revolute joint without limit specified
        else
            limit = [0,0]; % default fixed joint limit
        end
    end
    % Parent joint
    parent_joint_idx = idx_cell(child_link_names,parent_link_names{i_idx});
    if isempty(parent_joint_idx)
        parent_joint_name = '';
    else
        parent_joint_name = joint_names{parent_joint_idx};
    end
    % position and rotation offset
    if isfield(joint_i,'axis')
        a = cv(str2num(joint_i.axis.Attributes.xyz));
    else
        a = cv([0,0,0]);
    end
    % Joint position offset
    if isfield(joint_i,'origin')
        p_offset = cv(str2num(joint_i.origin.Attributes.xyz));
        R_offset = rpy2r(str2num(joint_i.origin.Attributes.rpy));
    else
        p_offset = cv([0,0,0]);
        R_offset = eye(3,3);
    end
    % Add joint
    q_rand = 0; dq = 0; ddq = 0;
    chain.joint(i_idx) = struct(...
        'name',joint_name,...
        'p',cv([0,0,0]),'R',eye(3,3),'a',cv(a),'type',joint_type,...
        'p_offset',cv(p_offset),'R_offset',R_offset,...
        'q',q_rand,'dq',dq,'ddq',ddq,'q_diff',0,'q_prev',0,...
        'v',cv([0,0,0]),'vo',cv([0,0,0]),'w',cv([0,0,0]),...
        'dvo',cv([0,0,0]),'dw',cv([0,0,0]),'u',0,'ext_f',cv([0,0,0]),...
        'parent',[],'childs',[],'link_idx',[],'limit',limit...
        );
    % Add parent joint index
    chain.joint(i_idx).parent = parent_joint_idx;
    fprintf('[%d/%d] joint name:[%s] type:[%s] added. \n',i_idx,n_joint,joint_name,joint_type);
end % for i_idx = 1:n_joint % for all joints

% Add the current joint as childs of the parent joint
for i_idx = 1:n_joint % for all joints
    joint_i = chain.joint(i_idx);
    parent_idx = joint_i.parent;
    if ~isempty(parent_idx)
        childs = chain.joint(parent_idx).childs;
        childs = [childs, i_idx]; % append
        chain.joint(parent_idx).childs = childs;
    end
end

%% Plot current chain
ca;
plot_chain(chain,'fig_pos',[0.0,0.3,0.4,0.65],...
    'PLOT_JOINT_SPHERE',1,'PLOT_JOINT_AXIS',1,'PLOT_JOINT_NAME',1);
plot_chain_graph(chain,'fig_pos',[0.0,0.0,0.4,0.3]);

%% Forward kinematics then plot
ca;
chain = fk_chain(chain,'');
plot_chain(chain,'fig_pos',[0.0,0.3,0.4,0.65],...
    'PLOT_JOINT_SPHERE',1,'PLOT_JOINT_AXIS',1,'PLOT_JOINT_NAME',1);

%% Parse link information
ca
n_link = length(robot.link); % number of link
chain.n_link = n_link;
chain.link_names = cell(1,n_link);
for i_idx = 1:chain.n_link % for all links
    link_i = robot.link{i_idx};
    link_name = link_i.Attributes.name;
    chain.link_names{i_idx} = link_name; % append link names
    fprintf('[%d/%d] link_name:[%s]\n',i_idx,chain.n_link,link_name);
end % for i_idx = 1:chain.n_link % for all links
fprintf('\n');
for i_idx = 1:chain.n_link % for all links
    link_i = robot.link{i_idx};
    link_name = link_i.Attributes.name;
    joint_idx = idx_cell(child_link_names,link_name); % attached joint index
    % Parse link mesh offset
    p_offset = cv([0,0,0]); R_offset = eye(3,3);
    if isfield(link_i,'visual')
        if isfield(link_i.visual,'origin')
            xyz = cv(str2num(link_i.visual.origin.Attributes.xyz));
            rpy = cv(str2num(link_i.visual.origin.Attributes.rpy));
            p_offset = xyz;
            R_offset = rpy2r(rpy);
        end
    end
    % Parse link mesh name and scale (if exist)
    filename = ''; mesh_path = ''; scale = cv([1,1,1]);
    if isfield(link_i,'visual')
        if isfield(link_i.visual,'geometry')
            if isfield(link_i.visual.geometry,'mesh')
                if isfield(link_i.visual.geometry.mesh.Attributes,'filename')
                    filename = link_i.visual.geometry.mesh.Attributes.filename;
                end
                if isfield(link_i.visual.geometry.mesh.Attributes,'scale')
                    scale = cv(str2num(link_i.visual.geometry.mesh.Attributes.scale));
                end
            end
        end
    end
    if ~isempty(filename) % if file name exists
        [~,name,ext] = fileparts(filename);
        [f,~,~] = fileparts(urdf_path);
        mesh_path = [f,'/visual/',name,ext];
    end
    if (~isempty(mesh_path)) && (~exist(mesh_path,'file'))
        % Check whether the specified stl file exist or not
        fprintf(2,'[%s] does not exist.\n',mesh_path);
    end
    % Load stl file
    fv = '';
    if exist(mesh_path,'file') % if file exists
        [~,~,ext] = fileparts(mesh_path);
        switch lower(ext)
            case '.stl'
                fv = load_stl(mesh_path);                 % load stl
                [fv.vertices,fv.faces]= patchslim(fv.vertices,fv.faces);    % reduce mesh
                fv.vertices = rv(scale).*fv.vertices;        % scaling
                fv.vertices = fv.vertices * R_offset';    % rotate mesh
                fv.vertices = fv.vertices + p_offset';    % translate mesh
            otherwise
                fprintf(2,'Unsupported file type:[%s]. Only .stl is supported. \n',ext);
        end
    end % if exist(chain.link(i_idx).mesh_path,'file') % if file exists
    
    % Parse box information
    box_size = ''; box_scale = cv([1,1,1]); box = '';
    if isfield(link_i,'visual')
        if isfield(link_i.visual,'geometry')
            if isfield(link_i.visual.geometry,'box')
                if isfield(link_i.visual.geometry.box.Attributes,'size')
                    box_size = cv(str2num(link_i.visual.geometry.box.Attributes.size));
                end
                if isfield(link_i.visual.geometry.box.Attributes,'scale')
                    box_scale = cv(str2num(link_i.visual.geometry.box.Attributes.scale));
                end
                if ~isempty(box_size)
                    box = struct('size',box_size,'scale',box_scale);
                end
            end
        end
    end
    % Append chain link
    chain.link(i_idx).name      = link_name;
    chain.link(i_idx).joint_idx = joint_idx;
    chain.link(i_idx).p_offset  = p_offset;
    chain.link(i_idx).R_offset  = R_offset;
    chain.link(i_idx).mesh_path = mesh_path;
    chain.link(i_idx).scale     = scale;
    chain.link(i_idx).fv        = fv; % append mesh
    chain.link(i_idx).box       = box;
    
    % Print
    fprintf('[%d/%d] link:[%s] added to joint:[%s]. mesh_path:[%s]\n',...
        i_idx,chain.n_link,link_name,chain.joint_names{joint_idx},mesh_path);
end % for i_idx = 1:chain.n_link % for all links

%% Plot chain with links added
ca;
plot_chain(chain,'fig_pos',[0.0,0.3,0.4,0.65],...
    'PLOT_JOINT_SPHERE',1,'PLOT_JOINT_AXIS',1,'PLOT_JOINT_NAME',1);

%% Optimize capsules
ca
SKIP_CAPSULE = 0;
for i_idx = 1:chain.n_link % for each link
    link_i = chain.link(i_idx);
    fv_i = link_i.fv;
    if (~isempty(fv_i)) && (~SKIP_CAPSULE)
        cap_opt_i = optimize_capsule(fv_i); % optimize capsule (takes time)
    else
        cap_opt_i = '';
    end
    chain.link(i_idx).capsule = cap_opt_i;
    fprintf('[%d/%d] capsule optimized\n',i_idx,chain.n_link);
end % for i_idx = 1:chain.n_link % for each link
fprintf('Done.\n');

%% Plot chain with capsules
ca;
plot_chain(chain,'fig_pos',[0.0,0.3,0.4,0.65],...
    'PLOT_JOINT_SPHERE',1,'PLOT_JOINT_AXIS',1,'PLOT_JOINT_NAME',1,...
    'PLOT_CAPSULE',1,'cfc','','cfa',0.5);

%% Update other information and plot
ca;
% Update mass, inertia, and com of links
chain = update_chain_mass_inertia_com(chain,'density',400);

% Forward kinematics
chain = fk_chain(chain,'');

% Get robot size
chain.sz = get_chain_sz(chain);
chain.axis_info = max(chain.sz.xyz_len)*[-1,+1,-1,+1,0,1.5]; % axis info

%%clc; ca;
tick = 0;
while true % loop
    tick = tick + 1;
    fig = plot_chain(chain,'fig_idx',1,'subfig_idx',1,'fig_pos',[0.5,0.5,0.4,0.5],...
        'view_info',[88,16],'axis_info',chain.axis_info,'AXIS_EQUAL',1,'USE_ZOOMRATE',1,...
        'PLOT_LINK',1,'llc','k','llw',2,'lls','-',...
        'PLOT_MESH',1,'mfc','','mfa',0.3,...
        'PLOT_BOX',1,'bfc','b','bfa',0.7,...
        'PLOT_CAPSULE',1,'cfc','','cfa',0.2,'cec','none',...
        'PLOT_BOX_ADDED',1,'bafa',0.4,...
        'PLOT_COM',1,'csc','r','csr',0.03,...
        'PLOT_JOINT_AXIS',1,'jal',0.05,'jalw',2,...
        'PLOT_JOINT_SPHERE',1,'jsr',0.02,'jsfc','k','jsfa',0.5,...
        'PLOT_ROTATE_AXIS',1,'ral',0.1,'raa',0.5,...
        'PLOT_JOINT_NAME',1,'jnfs',7 ...
        );
    plot_title(sprintf('tick:[%d]',tick),'fig_idx',1,'tfs',20,'interpreter','latex');
    drawnow; if ~ishandle(fig), break; end
end % while true % loop

%% Check self-collision
clc;ca;

% Get link indices to check self-collision
chain.sc_checks = get_sc_checks(chain,'collision_margin',max(chain.sz.xyz_len)/100);

% Check self-collision
while 1
    % Random robot pose
    q_rand = zeros(1,chain.n_rev_joint);
    for i_idx = 1:chain.n_rev_joint
        limit_i = chain.joint(idx_cell(chain.joint_names,chain.rev_joint_names{i_idx})).limit;
        q_rand(i_idx) = limit_i(1) + (limit_i(2)-limit_i(1))*rand;
    end
    chain = update_chain_q(chain,chain.rev_joint_names,q_rand);
    [SC,sc_link_pairs] = check_sc(chain,'collision_margin',0);
    if SC
        break;
    end
end

% Plot
plot_chain(chain,'fig_idx',1,'subfig_idx',1,'fig_pos',[0.0,0.5,0.4,0.5],...
    'PLOT_BOX',1,'bec','k');
for sc_idx = 1:size(sc_link_pairs,1)
    sc_link_pair = sc_link_pairs(sc_idx,:);
    i_idx = sc_link_pair(1); j_idx = sc_link_pair(2);
    % Link i
    [cap_i,T_i] = get_chain_capsule(chain,i_idx);
    % Link j
    [cap_j,T_j] = get_chain_capsule(chain,j_idx);
    % Plot colliding capsules
    plot_capsule(cap_i,'fig_idx',1,'subfig_idx',2*sc_idx-1,...
        'T',T_i,'cfc','r','cfa',0.5,'cec','k','cea',0.5);
    plot_capsule(cap_j,'fig_idx',1,'subfig_idx',2*sc_idx,...
        'T',T_j,'cfc','r','cfa',0.5,'cec','k','cea',0.5);
end
if SC
    title_str = sprintf('Self Collision Occurred');
    plot_title(title_str,'tfc','r','tfs',20);
else
    title_str = sprintf('No Collision');
    plot_title(title_str,'tfc','k','tfs',20);
end
fprintf('Done.\n');

%%



