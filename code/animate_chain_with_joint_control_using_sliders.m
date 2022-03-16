function animate_chain_with_joint_control_using_sliders(chain,varargin)
%
% Animate chain with joint control using sliders
%

% Parse input arguments
ps = inputParser;
addParameter(ps,'fig_pos_robot',[0.5,0.4,0.4,0.5]);
addParameter(ps,'fig_pos_slider',[0.5,0.2,0.4,0.2]);
addParameter(ps,'q_rev',get_q_rev_for_ik(chain));
addParameter(ps,'PLOT_COM',0);
addParameter(ps,'com_radius',0.05);
addParameter(ps,'PRINT_JOINT_POS',1);

addParameter(ps,'PLOT_MESH',true);                  % plot mesh
addParameter(ps,'mfc','');
addParameter(ps,'mfa',0.5);

addParameter(ps,'bfc','b');                         % box face color
addParameter(ps,'bafc','');                         % box added face color

addParameter(ps,'PLOT_LINK',true);                  % plot link
addParameter(ps,'PLOT_CAPSULE',true);               % plot capsule
addParameter(ps,'cfc',0.5*[1,1,1]);                 % capsule face color
addParameter(ps,'cfa',0.2);                         % capsule face alpha
addParameter(ps,'cec','none');                      % capsule edge color
addParameter(ps,'cea',0.5);                         % capsule edge alpha

addParameter(ps,'PLOT_JOINT_SPHERE',1);             %
addParameter(ps,'jsr',0.02);                        %

addParameter(ps,'PLOT_JOINT_AXIS',false);           % joint joint coordinates
addParameter(ps,'jal',0.1);                         % joint axis length
addParameter(ps,'jalw',2);                          % joint axis line width
addParameter(ps,'jals','-');                        % joint axis line style

addParameter(ps,'PLOT_ROTATE_AXIS',true);           % plot rotate axis
addParameter(ps,'ral',0.1);                         % rotate axis length
addParameter(ps,'rac','');                          % rotate axis color
addParameter(ps,'raa',0.5);                         % rotate axis alpha

addParameter(ps,'PLOT_JOINT_NAME',false);           % plot joint name
addParameter(ps,'PLOT_JOINT_TORQUE',false);         % plot joint torque
addParameter(ps,'jnfs',12);                         % joint name font size
addParameter(ps,'jnfn','consolas');                 % joint name font name
addParameter(ps,'jnfc','k');                        % joint name font color

addParameter(ps,'AXIS_OFF',0);
addParameter(ps,'NO_MARGIN',0);

addParameter(ps,'PLOT_GRAPH',0);

parse(ps,varargin{:});

fig_pos_robot       = ps.Results.fig_pos_robot;
fig_pos_slider      = ps.Results.fig_pos_slider;

PLOT_JOINT_NAME     = ps.Results.PLOT_JOINT_NAME;
PLOT_JOINT_TORQUE   = ps.Results.PLOT_JOINT_TORQUE;
jnfs                = ps.Results.jnfs;
jnfn                = ps.Results.jnfn;
jnfc                = ps.Results.jnfc;

q_rev               = ps.Results.q_rev;
PLOT_COM            = ps.Results.PLOT_COM;
com_radius          = ps.Results.com_radius;
PRINT_JOINT_POS     = ps.Results.PRINT_JOINT_POS;

PLOT_MESH           = ps.Results.PLOT_MESH;
mfc                 = ps.Results.mfc;
mfa                 = ps.Results.mfa;

bfc                 = ps.Results.bfc;
bafc                = ps.Results.bafc;

PLOT_LINK           = ps.Results.PLOT_LINK;
PLOT_CAPSULE        = ps.Results.PLOT_CAPSULE;
cfc                 = ps.Results.cfc;
cfa                 = ps.Results.cfa;
cec                 = ps.Results.cec;
cea                 = ps.Results.cea;

PLOT_JOINT_SPHERE   = ps.Results.PLOT_JOINT_SPHERE;
jsr                 = ps.Results.jsr;

PLOT_JOINT_AXIS     = ps.Results.PLOT_JOINT_AXIS;
jal                 = ps.Results.jal;
jalw                = ps.Results.jalw;
jals                = ps.Results.jals;

PLOT_ROTATE_AXIS    = ps.Results.PLOT_ROTATE_AXIS;
ral                 = ps.Results.ral;
rac                 = ps.Results.rac;
raa                 = ps.Results.raa;

AXIS_OFF            = ps.Results.AXIS_OFF;
NO_MARGIN           = ps.Results.NO_MARGIN;

PLOT_GRAPH          = ps.Results.PLOT_GRAPH;

% Constants
global g_cb         % global callback
global g_slider     % slider
R2D = 180/pi;
D2R = pi/180;

% Plot chain graph
if PLOT_GRAPH
    plot_chain_graph(chain,'fig_idx',3,'fig_pos',[0.5,0.5,0.4,0.3],...
        'interpreter','latex','text_fs',15,'max_str_len',inf,'title_str','','title_fs',30,...
        'NO_MARGIN',1,'x_margin',0.05,'y_margin',0.05);
end

% Figure with sliders
fig_slider = set_fig(figure(5),'pos',fig_pos_slider,...
    'HOLD_ON',0,'AXIS_EQUAL',0,'axis_info','',...
    'AXIS_OFF',0,'GRID_ON',0,'REMOVE_MENUBAR',1,'USE_DRAGZOOM',0,...
    'SET_CAMLIGHT',0,'SET_MATERIAL','','SET_AXISLABEL',0);
slider_names = chain.rev_joint_names;
slider_values = get_q_chain(chain,chain.rev_joint_names)*R2D; % to degree

bg_colors = zeros(chain.n_rev_joint,3);
for i_idx = 1:chain.n_rev_joint
    rev_joint_name = chain.rev_joint_names{i_idx};
    joint_idx = idx_cell(chain.joint_names,rev_joint_name);
    [~,max_idx] = max(abs(chain.joint(joint_idx).a));
    switch max_idx
        case 1
            bg_colors(i_idx,:) = [1,0.4,0.4];
        case 2
            bg_colors(i_idx,:) = [0.4,1,0.4];
        case 3
            bg_colors(i_idx,:) = [0.4,0.4,1];
    end
end
set_sliders_with_text_boxes(fig_slider,slider_names,slider_values,'bg_colors',bg_colors);

% Update slider val
for i_idx = 1:length(g_cb.slider_val)
    g_cb.slider_val{i_idx} = q_rev(i_idx)*R2D;
end

% Update chain self-collision checks
chain.sc_checks = get_sc_checks(chain);

% Animate robot model
tick = 0;
while true % loop
    % Update
    tick = tick + 1;
    for i_idx = 1:length(g_cb.slider_val)
        q_rev(i_idx) = g_cb.slider_val{i_idx}*D2R;
    end
    chain = update_chain_q(chain,chain.rev_joint_names,q_rev);
    SC = check_sc(chain);
    q_rev = get_q_chain(chain,chain.rev_joint_names); % get robot position
    for i_idx = 1:length(g_slider)
        if ishandle(g_slider{i_idx})
            set(g_slider{i_idx},'Value',q_rev(i_idx)*R2D); % apply limits to slider
        end
    end
    % Animate
    fig = plot_chain(chain,'fig_idx',1,'subfig_idx',1,'fig_pos',fig_pos_robot,...
        'view_info',[88,16],'axis_info',chain.axis_info,'AXIS_OFF',AXIS_OFF,'AXIS_EQUAL',1,'USE_ZOOMRATE',1,...
        'PLOT_LINK',PLOT_LINK,'llc','k','llw',1,'lls','-',...
        'PLOT_MESH',PLOT_MESH,'mfc',mfc,'mfa',mfa,...
        'PLOT_BOX',1,'bfc',bfc,'bfa',0.7,...
        'PLOT_CAPSULE',PLOT_CAPSULE,'cfc',cfc,'cfa',cfa,'cec',cec,'cea',cea,...
        'PLOT_BOX_ADDED',1,'bafc',bafc,'bafa',0.4,...
        'PLOT_COM',0,'csc','r','csr',0.01,...
        'PLOT_JOINT_AXIS',PLOT_JOINT_AXIS,'jal',jal,'jalw',jalw,'jals',jals,...
        'PLOT_JOINT_SPHERE',PLOT_JOINT_SPHERE,'jsr',jsr,'jsfc','k','jsfa',0.5,...
        'PLOT_ROTATE_AXIS',PLOT_ROTATE_AXIS,'ral',ral,'rac',rac,'raa',raa,...
        'PLOT_JOINT_NAME',PLOT_JOINT_NAME,'jnfs',jnfs,...
        'NO_MARGIN',NO_MARGIN...
        );
    if PLOT_COM
        com = get_chain_com(chain);
        plot_T(p2t(com),'fig_idx',1,'subfig_idx',1,'PLOT_AXIS',0,...
            'PLOT_SPHERE',1,'sfc','r','sr',com_radius,'sfa',0.9);
    end
    title_str = sprintf('tick:[%d] SC:[%d]',tick,SC);
    if SC
        tfc = 'r';
    else
        tfc = 'b';
    end
    plot_title_with_text(title_str,...
        'fig_idx',1,'tfs',20,'tfc',tfc,'interpreter','latex');
    drawnow; if ~ishandle(fig), break; end
    % Print-out joint position
    if PRINT_JOINT_POS
        fprintf('Rev. Joint: %s\n',...
            vec2str(rv(q_rev*R2D),'%+.1f'));
    end
end % while true % loop
if ~ishandle(fig), ca; end
fprintf('[animate_chain_with_joint_control_using_sliders] Done.\n');
