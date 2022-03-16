function chain = add_link_to_chain(chain,varargin)
%
% Add a link to a kinematic chain
%

% Parse options
ps = inputParser;
addParameter(ps,'name','');                     % link name
addParameter(ps,'joint_name','');               % joint name to attach link to
addParameter(ps,'p_offset',cv([0,0,0]));        %
addParameter(ps,'R_offset',eye(3,3));           %
addParameter(ps,'mesh_path','');                %
addParameter(ps,'scale',cv([1,1,1]));           %
addParameter(ps,'fv','');                       %
addParameter(ps,'box','');                      % box
addParameter(ps,'box_added','');                % additional box
addParameter(ps,'capsule','');                  % capsule
addParameter(ps,'v',cv([0,0,0]));               %
addParameter(ps,'vo',cv([0,0,0]));              %
addParameter(ps,'w',cv([0,0,0]));               %
parse(ps,varargin{:});
name                = ps.Results.name;
joint_name          = ps.Results.joint_name;
p_offset            = ps.Results.p_offset;
R_offset            = ps.Results.R_offset;
mesh_path           = ps.Results.mesh_path;
scale               = ps.Results.scale;
fv                  = ps.Results.fv;
box                 = ps.Results.box;
box_added           = ps.Results.box_added;
capsule             = ps.Results.capsule;
v                   = ps.Results.v;
vo                  = ps.Results.vo;
w                   = ps.Results.w;

% Parent joint index
joint_idx = idx_cell(chain.joint_names,joint_name);
if isempty(joint_idx)
    if ~isempty(joint_name)
        fprintf(2,'[add_link_to_chain] joint_name:[%s] is not in the chain.\n',...
            joint_name);
    end
end

% Add link
if ~isfield(chain,'link')
    % If link field does not exist
    chain.link(1) = struct( ...
        'name',name,'joint_idx',joint_idx, ...
        'p_offset',p_offset,'R_offset',R_offset,'mesh_path',mesh_path,'scale',scale,...
        'fv',fv,'box',box,'capsule',capsule,'box_added',box_added, ...
        'm',0,'I_bar',eye(3,3),'com_bar',cv([0,0,0]),...
        'v',v,'vo',vo,'w',w ...
        );
    chain.link_names{1} = name;
    chain.n_link = 1;
else
    % Add link
    chain.n_link = chain.n_link + 1;
    chain.link_names{chain.n_link} = name;
    chain.link(chain.n_link) = struct( ...
        'name',name,'joint_idx',joint_idx, ...
        'p_offset',p_offset,'R_offset',R_offset,'mesh_path',mesh_path,'scale',scale,...
        'fv',fv,'box',box,'capsule',capsule,'box_added',box_added, ...
        'm',0,'I_bar',eye(3,3),'com_bar',cv([0,0,0]),...
        'v',v,'vo',vo,'w',w ...
        );
end

% Update joint information (there could be only one link index per joint)
if ~isempty(joint_idx)
    chain.joint(joint_idx).link_idx = chain.n_link;
end
