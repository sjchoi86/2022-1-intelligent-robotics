function chain = get_chain_from_urdf_with_caching(robot_name,varargin)
%
% Get a kinematic chain from URDF
%

% Parse options
ps = inputParser;
addParameter(ps,'RE',false); % redo loading
addParameter(ps,'SKIP_CAPSULE',false);
addParameter(ps,'urdf_path',sprintf('../urdf/%s/%s_urdf.xml',robot_name,robot_name));
addParameter(ps,'cache_folder','../cache'); % folder to caache
parse(ps,varargin{:});
RE           = ps.Results.RE;
urdf_path    = ps.Results.urdf_path;
cache_folder = ps.Results.cache_folder;
SKIP_CAPSULE = ps.Results.SKIP_CAPSULE;

% Cache file (fixed file path)
cache_path = sprintf('%s/model/urdf_%s.mat',cache_folder,robot_name);
[p,~,~] = fileparts(cache_path);
make_dir_if_not_exist(p);

if exist(cache_path,'file') && (RE==0)
    l = load(cache_path); % load
    chain = l.chain;
else
    chain = get_chain_from_urdf(robot_name,...
        'urdf_path',urdf_path,'SKIP_CAPSULE',SKIP_CAPSULE);
    save(cache_path,'chain');
    fprintf(2,'[%s] saved.\n',cache_path);
end

% Post processing 
chain = update_chain_q(chain,chain.rev_joint_names,zeros(chain.n_rev_joint,1),...
    'IGNORE_LIMIT',0,'FK',1);
