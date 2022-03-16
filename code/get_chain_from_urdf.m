function chain = get_chain_from_urdf(robot_name,varargin)
%
% Get a kinematic chain from URDF
%

% Parse options
ps = inputParser;
addParameter(ps,'density',500);    % density (e.g., 500 [kg/m^3] for firm wood)
addParameter(ps,'urdf_path',sprintf('../urdf/%s/%s_urdf.xml',robot_name,robot_name));
addParameter(ps,'SKIP_CAPSULE',0);
parse(ps,varargin{:});
density         = ps.Results.density;
urdf_path       = ps.Results.urdf_path;
SKIP_CAPSULE    = ps.Results.SKIP_CAPSULE;

% Parse xml 
s = xml2str(urdf_path);
robot = s.robot;

% Parse information
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

% Parse joint information
for i_idx = 1:n_joint % for all joints
    joint_i = robot.joint{i_idx};
    joint_name = joint_i.Attributes.name;
    joint_type = joint_i.Attributes.type;
    % Rotation axis
    if isequal(joint_type,'revolute')
        a = cv(str2num(joint_i.axis.Attributes.xyz));
    else
        a = cv([0,0,0]);
    end
    % Joint position limit
    if isfield(joint_i,'limit')
        limit_lower = joint_i.limit.Attributes.lower;
        upper = joint_i.limit.Attributes.upper;
        if limit_lower(1) == '$'
            limit_lower = eval(limit_lower(3:end-1));
        else
            limit_lower = eval(limit_lower);
        end
        if upper(1) == '$'
            upper = eval(upper(3:end-1));
        else
            upper = eval(upper);
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
    q = 0; dq = 0; ddq = 0;
    chain.joint(i_idx) = struct(...
        'name',joint_name,...
        'p',cv([0,0,0]),'R',eye(3,3),'a',cv(a),'type',joint_type,...
        'p_offset',cv(p_offset),'R_offset',R_offset,...
        'q',q,'dq',dq,'ddq',ddq,'q_diff',0,'q_prev',0,...
        'v',cv([0,0,0]),'vo',cv([0,0,0]),'w',cv([0,0,0]),...
        'dvo',cv([0,0,0]),'dw',cv([0,0,0]),...
        'u',0,'ext_f',cv([0,0,0]),...
        'parent',[],'childs',[],...
        'link_idx',[],'limit',limit...
        );
    % Add parent joint index
    chain.joint(i_idx).parent = parent_joint_idx;
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

% Forward kinematics
chain = fk_chain(chain,'');

% Parse link information
n_link = length(robot.link); % number of link
for i_idx = 1:n_link % for all links
    link_i = robot.link{i_idx};
    link_name = link_i.Attributes.name;
    joint_idx = idx_cell(child_link_names,link_name); % attached joint index
    if isempty(joint_idx)
        joint_name = '';
    else
        joint_name = chain.joint_names{joint_idx};
    end
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
    
    % Optimize capsule
    if (~isempty(fv)) && (~SKIP_CAPSULE)
        cap_opt = optimize_capsule(fv); % optimize capsule (takes time)
    else
        cap_opt = '';
    end
    
    % Append chain link
    chain = add_link_to_chain(chain,'name',link_name,'joint_name',joint_name,...
        'p_offset',p_offset,'R_offset',R_offset,'mesh_path',mesh_path,'scale',scale,'fv',fv,...
        'box',box,'box_added','','capsule',cap_opt,...
        'v',cv([0,0,0]),'vo',cv([0,0,0]),'w',cv([0,0,0]));

end % for i_idx = 1:chain.n_link % for all links

% Update mass, inertia, and com of links
chain = update_chain_mass_inertia_com(chain,'density',density);

% Forward kinematics
chain = fk_chain(chain,'');

% Get robot size
chain.sz = get_chain_sz(chain);
chain.axis_info = max(chain.sz.xyz_len)*[-1,+1,-1,+1,0,1.5]; % axis info
