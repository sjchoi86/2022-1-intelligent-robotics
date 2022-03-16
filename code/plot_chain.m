function fig = plot_chain(chain,varargin)
%
% Plot a kinematic chain
%
% % Basic usage
%
% tick = 0;
% while 1
%     tick = tick + 1;
%     fig = plot_chain(chain,'fig_idx',1,'subfig_idx',1,'fig_pos',[0.5,0.5,0.4,0.5],...
%         'view_info',[88,16],'axis_info',chain.axis_info,'AXIS_EQUAL',1,'USE_ZOOMRATE',1,...
%         'PLOT_LINK',1,'llc','k','llw',2,'lls','-',...
%         'PLOT_MESH',1,'mfc','','mfa',0.3,...
%         'PLOT_BOX',1,'bfc','b','bfa',0.7,...
%         'PLOT_CAPSULE',1,'cfc','','cfa',0.2,'cec','none',...
%         'PLOT_BOX_ADDED',1,'bafa',0.4,...
%         'PLOT_COM',1,'csc','r','csr',0.03,...
%         'PLOT_JOINT_AXIS',1,'jal',0.05,'jalw',2,...
%         'PLOT_JOINT_SPHERE',1,'jsr',0.02,'jsfc','k','jsfa',0.5,...
%         'PLOT_ROTATE_AXIS',1,'ral',0.1,'raa',0.5,...
%         'PLOT_JOINT_NAME',1,'jnfs',7 ...
%         );
%     plot_title(sprintf('tick:[%d]',tick),'fig_idx',1,'tfs',20,'interpreter','latex');
%     drawnow; if ~ishandle(fig), break; end
% end
%

global zoomPct % zoom rate from dragzoom
persistent h

% Make enough handlers at the first
if isempty(h)
    for i = 1:10
        for j = 1:100
            h{i,j}.first_flag = true;
        end
    end
end

% Parse options
ps = inputParser;
addParameter(ps,'fig_idx',1);
addParameter(ps,'subfig_idx',1);
addParameter(ps,'fig_pos',[0.5,0.4,0.5,0.6]);       % position of the figure
addParameter(ps,'view_info',[88,16]);               % view information
addParameter(ps,'axis_info','');                    % axis information
addParameter(ps,'AXIS_EQUAL',1);                    % axis equal
addParameter(ps,'AXIS_OFF',0);                      % axis off
addParameter(ps,'GRID_ON',1);                       % grid on
addParameter(ps,'REMOVE_MENUBAR',1);                % remove menubar
addParameter(ps,'USE_DRAGZOOM',1);                  % use dragzoom
addParameter(ps,'USE_ZOOMRATE',1);                  % interactively change fontsize and linewidth
addParameter(ps,'gcf_color','w');                   % background color

addParameter(ps,'SET_CAMLIGHT',1);                  % set camera light
addParameter(ps,'SET_MATERIAL','METAL');            % set material property 'SHINY' 'DULL' 'METAL'
addParameter(ps,'SET_AXISLABEL',1);                 % axis label
addParameter(ps,'afs',13);                          % axis font size
addParameter(ps,'interpreter','latex');             % interpreter

addParameter(ps,'PLOT_LINK',true);                  % plot link
addParameter(ps,'llc','k');                         % link line color
addParameter(ps,'llw',2);                           % link line width
addParameter(ps,'lls','-');                         % link line style

addParameter(ps,'PLOT_MESH',true);                  % plot mesh
addParameter(ps,'mfc',0.5*[1,1,1]);                 % mesh face color
addParameter(ps,'mfa',0.4);                         % mesh face alpha

addParameter(ps,'PLOT_BOX',true);                   % plot box
addParameter(ps,'bfc','b');                         % box face color
addParameter(ps,'bfa',0.5);                         % box face alpha
addParameter(ps,'bec','none');                      % box edge color

addParameter(ps,'PLOT_CAPSULE',false);              % plot capsule
addParameter(ps,'cfc',[0.7,0.7,0.99]);              % capsule face color
addParameter(ps,'cfa',0.2);                         % capsule face alpha
addParameter(ps,'cec','none');                      % capsule edge color
addParameter(ps,'cea',0.5);                         % capsule edge alpha

addParameter(ps,'PLOT_BOX_ADDED',true);             % plot box added
addParameter(ps,'bafc','');                         % box added face color
addParameter(ps,'bafa','');                         % box added face alpha
addParameter(ps,'baec','');                         % box added edge color

addParameter(ps,'PLOT_COM',false);                  % plot com
addParameter(ps,'csc','r');                         % com sphere color
addParameter(ps,'csr',0.1);                         % com sphere radius
addParameter(ps,'csa',0.5);                         % com sphere alphs

addParameter(ps,'PLOT_LINK_V',false);               % plot com directional velocity
addParameter(ps,'lvfc','m');                        % com directional velocity face color
addParameter(ps,'lvfa',0.5);                        % com directional velocity face alpha
addParameter(ps,'lvar',0.05);                       % com directional velocity arrow rate
addParameter(ps,'lvsw',0.05);                       % com directional velocity stem width
addParameter(ps,'lvtw',0.1);                        % com directional velocity tip width

addParameter(ps,'PLOT_LINK_W',false);               % plot com angular velocity
addParameter(ps,'lwfc','c');                        % com angular velocity face color
addParameter(ps,'lwfa',0.5);                        % com angular velocity face alpha
addParameter(ps,'lwar',0.05);                       % com angular velocity arrow rate
addParameter(ps,'lwsw',0.05);                       % com angular velocity stem width
addParameter(ps,'lwtw',0.1);                        % com angular velocity tip width

addParameter(ps,'PLOT_JOINT_AXIS',false);           % joint joint coordinates
addParameter(ps,'jal',0.1);                         % joint axis length
addParameter(ps,'jalw',2);                          % joint axis line width
addParameter(ps,'jals','-');                        % joint axis line style

addParameter(ps,'PLOT_JOINT_SPHERE',false);         % plot joint with a sphere
addParameter(ps,'jsr',0.02);                        % joint sphere radius
addParameter(ps,'jsfc','k');                        % joint sphere face color
addParameter(ps,'jsfa',0.5);                        % joint sphere face alpha

addParameter(ps,'PLOT_ROTATE_AXIS',true);           % plot rotate axis
addParameter(ps,'ral',0.1);                         % rotate axis length
addParameter(ps,'rac','');                          % rotate axis color
addParameter(ps,'raa',0.5);                         % rotate axis alpha

addParameter(ps,'PLOT_JOINT_NAME',false);           % plot joint name
addParameter(ps,'PLOT_JOINT_TORQUE',false);         % plot joint torque
addParameter(ps,'jnfs',12);                         % joint name font size
addParameter(ps,'jnfn','consolas');                 % joint name font name
addParameter(ps,'jnfc','k');                        % joint name font color

addParameter(ps,'PLOT_JOINT_V',false);              % plot joint directional velocity
addParameter(ps,'jvfc','m');                        % joint directional velocity face color
addParameter(ps,'jvfa',0.5);                        % joint directional velocity face alpha
addParameter(ps,'jvar',0.05);                       % joint directional velocity arrow rate
addParameter(ps,'jvsw',0.05);                       % joint directional velocity stem width
addParameter(ps,'jvtw',0.1);                        % joint directional velocity tip width

addParameter(ps,'PLOT_JOINT_W',false);              % plot joint angular velocity
addParameter(ps,'jwfc','c');                        % joint angular velocity face color
addParameter(ps,'jwfa',0.5);                        % joint angular velocity face alpha
addParameter(ps,'jwar',0.05);                       % joint angular velocity arrow rate
addParameter(ps,'jwsw',0.05);                       % joint angular velocity stem width
addParameter(ps,'jwtw',0.1);                        % joint angular velocity tip width

addParameter(ps,'NO_MARGIN',0);                     % no margin on axes

parse(ps,varargin{:});

fig_idx             = ps.Results.fig_idx;
subfig_idx          = ps.Results.subfig_idx;
fig_pos             = ps.Results.fig_pos;
view_info           = ps.Results.view_info;
axis_info           = ps.Results.axis_info;
AXIS_EQUAL          = ps.Results.AXIS_EQUAL;
AXIS_OFF            = ps.Results.AXIS_OFF;
GRID_ON             = ps.Results.GRID_ON;
REMOVE_MENUBAR      = ps.Results.REMOVE_MENUBAR;
USE_DRAGZOOM        = ps.Results.USE_DRAGZOOM;
USE_ZOOMRATE        = ps.Results.USE_ZOOMRATE;
gcf_color           = ps.Results.gcf_color;

SET_CAMLIGHT        = ps.Results.SET_CAMLIGHT;
SET_MATERIAL        = ps.Results.SET_MATERIAL;
SET_AXISLABEL       = ps.Results.SET_AXISLABEL;
afs                 = ps.Results.afs;
interpreter         = ps.Results.interpreter;

PLOT_LINK           = ps.Results.PLOT_LINK;
llc                 = ps.Results.llc;
llw                 = ps.Results.llw;
lls                 = ps.Results.lls;

PLOT_BOX            = ps.Results.PLOT_BOX;
bfc                 = ps.Results.bfc;
bfa                 = ps.Results.bfa;
bec                 = ps.Results.bec;

PLOT_MESH           = ps.Results.PLOT_MESH;
mfc                 = ps.Results.mfc;
mfa                 = ps.Results.mfa;

PLOT_CAPSULE        = ps.Results.PLOT_CAPSULE;
cfc                 = ps.Results.cfc;
cfa                 = ps.Results.cfa;
cec                 = ps.Results.cec;
cea                 = ps.Results.cea;

PLOT_BOX_ADDED      = ps.Results.PLOT_BOX_ADDED;
bafc                = ps.Results.bafc;
bafa                = ps.Results.bafa;
baec                = ps.Results.baec;

PLOT_COM            = ps.Results.PLOT_COM;
csc                 = ps.Results.csc;
csr                 = ps.Results.csr;
csa                 = ps.Results.csa;

PLOT_LINK_V         = ps.Results.PLOT_LINK_V;
lvfc                = ps.Results.lvfc;
lvfa                = ps.Results.lvfa;
lvar                = ps.Results.lvar;
lvsw                = ps.Results.lvsw;
lvtw                = ps.Results.lvtw;

PLOT_LINK_W         = ps.Results.PLOT_LINK_W;
lwfc                = ps.Results.lwfc;
lwfa                = ps.Results.lwfa;
lwar                = ps.Results.lwar;
lwsw                = ps.Results.lwsw;
lwtw                = ps.Results.lwtw;

PLOT_JOINT_AXIS     = ps.Results.PLOT_JOINT_AXIS;
jal                 = ps.Results.jal;
jalw                = ps.Results.jalw;
jals                = ps.Results.jals;

PLOT_JOINT_SPHERE   = ps.Results.PLOT_JOINT_SPHERE;
jsr                 = ps.Results.jsr;
jsfc                = ps.Results.jsfc;
jsfa                = ps.Results.jsfa;

PLOT_ROTATE_AXIS    = ps.Results.PLOT_ROTATE_AXIS;
ral                 = ps.Results.ral;
rac                 = ps.Results.rac;
raa                 = ps.Results.raa;
rasw                = ral/8;
ratw                = ral/4;

PLOT_JOINT_NAME     = ps.Results.PLOT_JOINT_NAME;
PLOT_JOINT_TORQUE   = ps.Results.PLOT_JOINT_TORQUE;
jnfs                = ps.Results.jnfs;
jnfn                = ps.Results.jnfn;
jnfc                = ps.Results.jnfc;

PLOT_JOINT_V        = ps.Results.PLOT_JOINT_V;
jvfc                = ps.Results.jvfc;
jvfa                = ps.Results.jvfa;
jvar                = ps.Results.jvar;
jvsw                = ps.Results.jvsw;
jvtw                = ps.Results.jvtw;

PLOT_JOINT_W        = ps.Results.PLOT_JOINT_W;
jwfc                = ps.Results.jwfc;
jwfa                = ps.Results.jwfa;
jwar                = ps.Results.jwar;
jwsw                = ps.Results.jwsw;
jwtw                = ps.Results.jwtw;

NO_MARGIN           = ps.Results.NO_MARGIN;

% Dragzoom rate
if isempty(zoomPct) || (USE_ZOOMRATE==0)
    dragzoom_rate = 1;
else
    dragzoom_rate = 1 + (zoomPct-100)/150; % normalize by 100
    max_dragzoom_rate = 8;
    if dragzoom_rate > max_dragzoom_rate
        dragzoom_rate = max_dragzoom_rate;
    elseif dragzoom_rate < 1/max_dragzoom_rate
        dragzoom_rate = 1/max_dragzoom_rate;
    end
end

if h{fig_idx,subfig_idx}.first_flag || ...
        (~ishandle(h{fig_idx,subfig_idx}.fig))
    % First time
    h{fig_idx,subfig_idx}.first_flag = false;
    
    % Set figure
    h{fig_idx,subfig_idx}.fig = set_fig(figure(fig_idx),...
        'pos',fig_pos,'view_info',view_info,'axis_info',axis_info,'AXIS_EQUAL',AXIS_EQUAL,...
        'GRID_ON',GRID_ON,'REMOVE_MENUBAR',REMOVE_MENUBAR,'USE_DRAGZOOM',USE_DRAGZOOM,...
        'SET_CAMLIGHT',SET_CAMLIGHT,'SET_MATERIAL',SET_MATERIAL,'SET_AXISLABEL',SET_AXISLABEL,'AXIS_OFF',AXIS_OFF,...
        'afs',afs,'interpreter',interpreter,'gcf_color',gcf_color,'NO_MARGIN',NO_MARGIN ...
        );
    
    % Plot link
    if PLOT_LINK
        for i_idx = 1:chain.n_joint
            joint_i = chain.joint(i_idx);
            parent = joint_i.parent;
            if ~isempty(parent)
                joint_fr = chain.joint(parent);
                p1 = joint_fr.p;
                p2 = joint_i.p;
                h{fig_idx,subfig_idx}.line{i_idx} = plot3(...
                    [p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],...
                    'Color',llc,'LineWidth',llw,'LineStyle',lls);
            end
        end
    end
    
    % Plot mesh
    if PLOT_MESH && isfield(chain,'link')
        colors = linspecer(chain.n_link);
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if isfield(link_i,'fv') % if 'fv' exists
                joint_idx = link_i.joint_idx;
                if ~isempty(joint_idx)
                    % Mother joint position
                    p_link = chain.joint(joint_idx).p;
                    R_link = chain.joint(joint_idx).R;
                    % Plot mesh
                    fv = link_i.fv;
                    if (~isempty(fv)) % if mesh exists
                        % Plot mesh with patch
                        V = fv.vertices;
                        V = V*(1.0 + 0.01*subfig_idx); % little size trick
                        if isempty(mfc)
                            mfc_i = colors(i_idx,:);
                        else
                            mfc_i = mfc;
                        end
                        h{fig_idx,subfig_idx}.mesh{i_idx} = ...
                            patch('faces',fv.faces,'vertices',V,...
                            'FaceColor', mfc_i, ...
                            'EdgeColor', 'none','FaceLighting','gouraud',...
                            'AmbientStrength', 0.2, 'FaceAlpha', mfa);
                        % Define hgtransform
                        h{fig_idx,subfig_idx}.mesh_t{i_idx} = hgtransform;
                        set(h{fig_idx,subfig_idx}.mesh{i_idx},...
                            'parent',h{fig_idx,subfig_idx}.mesh_t{i_idx});
                        % Move
                        tform = pr2t(p_link,R_link);
                        set(h{fig_idx,subfig_idx}.mesh_t{i_idx},'Matrix',tform);
                    end
                end
            end % if isfield(link_i,'fv') % if 'fv' exists
        end % for i_idx = 1:chain.n_link % for all links
    end % if PLOT_MESH && isfield(chain,'link')
    
    % Plot box
    if PLOT_BOX && isfield(chain,'link')
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            box_i = link_i.box;
            if ~isempty(box_i) && ~isempty(link_i.joint_idx) % if box exists
                box_size  = box_i.size;
                box_scale = box_i.scale;
                joint_idx = link_i.joint_idx;
                p_link = chain.joint(joint_idx).p;
                R_link = chain.joint(joint_idx).R;
                p_offset = link_i.p_offset;
                R_offset = link_i.R_offset;
                xyz_min = p_offset-box_size.*box_scale/2; % box at the center
                xyz_len = box_size.*box_scale;
                vertex_matrix = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
                faces_matrix = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
                vertex_matrix = vertex_matrix.*xyz_len';
                vertex_matrix = vertex_matrix + xyz_min'; % do basic translation
                vertex_matrix = vertex_matrix * R_offset';
                h{fig_idx,subfig_idx}.box{i_idx} = ...
                    patch('Vertices',vertex_matrix,'Faces',faces_matrix,...
                    'FaceColor',bfc,'FaceAlpha',bfa,'EdgeColor',bec,...
                    'facelighting','gouraud');
                h{fig_idx,subfig_idx}.box_t{i_idx} = hgtransform;
                set(h{fig_idx,subfig_idx}.box{i_idx},...
                    'parent',h{fig_idx,subfig_idx}.box_t{i_idx});
                tform = pr2t(p_link,R_link);
                set(h{fig_idx,subfig_idx}.box_t{i_idx},'Matrix',tform);
            end
        end % for i_idx = 1:chain.n_link % for all links
    end
    
    % Plot box added
    if PLOT_BOX_ADDED && isfield(chain,'link')
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if isfield(link_i,'box_added') % if 'box_added' exists
                if ~isempty(link_i.box_added)
                    box_added = link_i.box_added;
                    xyz_min = cv(box_added.xyz_min);
                    xyz_len = cv(box_added.xyz_len);
                    vertex_matrix = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
                    faces_matrix = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
                    vertex_matrix = vertex_matrix.*xyz_len';
                    vertex_matrix = vertex_matrix + xyz_min';
                    if ~isfield(box_added,'ec')
                        box_added.ec = 'k';
                    end
                    if ~isempty(baec)
                        box_added.ec = baec;
                    end
                    if isempty(bafa)
                        bafa_use = box_added.alpha;
                    else
                        bafa_use = bafa;
                    end
                    if isempty(bafc)
                        bafc_use = box_added.color;
                    else
                        bafc_use = bafc;
                    end
                    h{fig_idx,subfig_idx}.box_added{i_idx} = patch(...
                        'Vertices',vertex_matrix,...
                        'Faces',faces_matrix,...
                        'FaceColor',bafc_use,'FaceAlpha',bafa_use,...
                        'EdgeColor',box_added.ec,'lineWidth',1);
                    h{fig_idx,subfig_idx}.box_added_t{i_idx} = hgtransform;
                    set(h{fig_idx,subfig_idx}.box_added{i_idx},...
                        'parent',h{fig_idx,subfig_idx}.box_added_t{i_idx});
                    T = pr2t(chain.joint(link_i.joint_idx).p, chain.joint(link_i.joint_idx).R);
                    if isfield(box_added,'p_offset')
                        p_offset = box_added.p_offset;
                        R_offset = box_added.R_offset;
                        T_offset = pr2t(p_offset,R_offset);
                    else
                        T_offset = pr2t(cv([0,0,0]),eye(3,3));
                    end
                    tform = T*T_offset;
                    set(h{fig_idx,subfig_idx}.box_added_t{i_idx},'Matrix',tform);
                end
            end % if isfield(link_i,'box_added') % if 'box_added' exists
        end % for i_idx = 1:chain.n_link % for all links
    end % if PLOT_BOX_ADDED && isfield(chain,'link')
    
    % Plot capsule
    if PLOT_CAPSULE && isfield(chain,'link')
        colors = linspecer(chain.n_link);
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.capsule)
                cap = link_i.capsule;
                joint_idx = link_i.joint_idx;
                if isempty(joint_idx)
                    joint_idx = get_topmost_idx(chain);
                end
                p_link = chain.joint(joint_idx).p;
                R_link = chain.joint(joint_idx).R;
                if isempty(cfc)
                    cfc_i = colors(i_idx,:);
                else
                    cfc_i = cfc;
                end
                if isempty(cec)
                    cec_i = colors(i_idx,:);
                else
                    cec_i = cec;
                end
                [p_offset,R_offset] = t2pr(cap.T_offset);
                vertices = rv(p_offset) + cap.vertices*R_offset';
                h{fig_idx,subfig_idx}.cap_patch{i_idx} = ...
                    patch('faces',cap.faces,'vertices',vertices,...
                    'FaceColor',cfc_i,'EdgeColor',cec_i,'EdgeAlpha',cea,...
                    'FaceLighting','gouraud','AmbientStrength', 0.5, 'FaceAlpha', cfa);
                h{fig_idx,subfig_idx}.cap_patch_t{i_idx} = hgtransform;
                set(h{fig_idx,subfig_idx}.cap_patch{i_idx},...
                    'parent',h{fig_idx,subfig_idx}.cap_patch_t{i_idx});
                tform = pr2t(p_link,R_link);
                
                set(h{fig_idx,subfig_idx}.cap_patch_t{i_idx},'Matrix',tform);
            end
        end % for i_idx = 1:chain.n_link % for all links
    end % if PLOT_BOX_ADDED && isfield(chain,'link')
    
    % Plot com
    if PLOT_COM && isfield(chain,'link')
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.com_bar) && ~isempty(link_i.joint_idx)
                joint_idx = link_i.joint_idx;
                com_bar = link_i.com_bar; % in local coordinates (post-multiply this)
                p_joint = chain.joint(joint_idx).p;
                R_joint = chain.joint(joint_idx).R;
                T_joint = pr2t(p_joint,R_joint);
                [x,y,z] = ellipsoid(0,0,0,csr,csr,csr,30);
                fv = surf2patch(x,y,z);
                h{fig_idx,subfig_idx}.com{i_idx} = patch(fv,...
                    'EdgeColor','none','FaceColor',csc,'FaceAlpha',csa,'facelighting','gouraud');
                h{fig_idx,subfig_idx}.com_t{i_idx} = hgtransform;
                set(h{fig_idx,subfig_idx}.com{i_idx},'parent',h{fig_idx,subfig_idx}.com_t{i_idx});
                tform = T_joint*pr2t(com_bar,''); % com pose
                set(h{fig_idx,subfig_idx}.com_t{i_idx},'Matrix',tform);
            end
        end % for i_idx = 1:chain.n_link % for all links
    end % if PLOT_COM && isfield(chain,'link')
    
    % Plot link directional velocity
    if PLOT_LINK_V
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.com_bar)
                joint_idx = link_i.joint_idx;
                com_bar = link_i.com_bar; % in local coordinates (post-multiply this)
                p_joint = chain.joint(joint_idx).p;
                R_joint = chain.joint(joint_idx).R;
                T_joint = pr2t(p_joint,R_joint);
                T_com = T_joint*pr2t(com_bar,'');
                p_com = t2p(T_com);
                lv_uv = uv(link_i.v);          % link directional velocity unit vector
                lv_len = jvar*norm(link_i.v);  % link directional velocity arrow length
                p1 = p_com;
                p2 = p_com + lv_len*lv_uv;
                fv = get_arrow_3d(p1,p2,'color',lvfc,'stemWidth',lvsw,'tipWidth',lvtw,...
                    'facealpha',lvfa);
                h{fig_idx,subfig_idx}.lv_arrow{i_idx} = patch(fv,...
                    'facecolor',lvfc,'edgeColor','none','FaceAlpha',lvfa,'FaceLighting','gouraud');
            end
        end
    end
    
    % Plot link angular velocity
    if PLOT_LINK_W
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.com_bar)
                joint_idx = link_i.joint_idx;
                com_bar = link_i.com_bar; % in local coordinates (post-multiply this)
                p_joint = chain.joint(joint_idx).p;
                R_joint = chain.joint(joint_idx).R;
                T_joint = pr2t(p_joint,R_joint);
                T_com = T_joint*pr2t(com_bar,'');
                p_com = t2p(T_com);
                lw_uv = uv(link_i.w);          % link directional velocity unit vector
                lw_len = jvar*norm(link_i.w);  % link directional velocity arrow length
                p1 = p_com;
                p2 = p_com + lw_len*lw_uv;
                fv = get_arrow_3d(p1,p2,'color',lwfc,'stemWidth',lwsw,'tipWidth',lwtw,...
                    'facealpha',lwfa);
                h{fig_idx,subfig_idx}.lw_arrow{i_idx} = patch(fv,...
                    'facecolor',lwfc,'edgeColor','none','FaceAlpha',lwfa,'FaceLighting','gouraud');
            end
        end
    end
    
    % Plot joint related information
    for i_idx = 1:chain.n_joint % plot joints
        joint_i = chain.joint(i_idx);
        p = joint_i.p;
        R = joint_i.R;
        ex = R(:,1); ey = R(:,2); ez = R(:,3);
        
        if PLOT_JOINT_AXIS
            h{fig_idx,subfig_idx}.cx{i_idx} = plot3([p(1),p(1)+jal*ex(1)],[p(2),p(2)+jal*ex(2)],...
                [p(3),p(3)+jal*ex(3)],'r','LineWidth',jalw,'LineStyle',jals);
            h{fig_idx,subfig_idx}.cy{i_idx} = plot3([p(1),p(1)+jal*ey(1)],[p(2),p(2)+jal*ey(2)],...
                [p(3),p(3)+jal*ey(3)],'g','LineWidth',jalw,'LineStyle',jals);
            h{fig_idx,subfig_idx}.cz{i_idx} = plot3([p(1),p(1)+jal*ez(1)],[p(2),p(2)+jal*ez(2)],...
                [p(3),p(3)+jal*ez(3)],'b','LineWidth',jalw,'LineStyle',jals);
        end
        
        if PLOT_ROTATE_AXIS
            a = joint_i.a; % axis
            if ~isequal(a,cv([0,0,0]))
                p1 = p; p2 = p+R*a*ral;
                if isempty(rac)
                    [~,max_idx] = max(abs(a'));
                    switch max_idx
                        case 1, rac_i = 'r';
                        case 2, rac_i = 'g';
                        case 3, rac_i = 'b';
                    end
                else
                    rac_i = rac;
                end
                fv = get_arrow_3d(p1,p2,'color',rac_i,'stemWidth',rasw,'tipWidth',ratw,...
                    'facealpha',raa);
                h{fig_idx,subfig_idx}.arrow{i_idx} = patch(fv,'facecolor',rac_i,'edgeColor','none',...
                    'FaceAlpha',raa,'FaceLighting','gouraud');
            end
        end
        
        if PLOT_JOINT_SPHERE
            [x,y,z] = ellipsoid(0,0,0,jsr,jsr,jsr,30);
            fv = surf2patch(x,y,z);
            h{fig_idx,subfig_idx}.sphere{i_idx} = patch(fv,...
                'EdgeColor','none','FaceColor',jsfc,'FaceAlpha',jsfa,'facelighting','gouraud');
            h{fig_idx,subfig_idx}.sphere_t{i_idx} = hgtransform;
            set(h{fig_idx,subfig_idx}.sphere{i_idx},'parent',h{fig_idx,subfig_idx}.sphere_t{i_idx});
            tform = pr2t(p,'');
            set(h{fig_idx,subfig_idx}.sphere_t{i_idx},'Matrix',tform);
        end
        
        if PLOT_JOINT_NAME || PLOT_JOINT_TORQUE
            text_p = p;
            if PLOT_JOINT_TORQUE && isequal(joint_i.type,'revolute')
                text_str = sprintf(' [%d]%s(%.2f)',i_idx,joint_i.name,joint_i.u);
            else
                text_str = sprintf(' [%d]%s',i_idx,joint_i.name);
            end
            h{fig_idx,subfig_idx}.text{i_idx} = text(text_p(1),text_p(2),text_p(3),...
                text_str,...
                'FontSize',jnfs,'FontName',jnfn,'Color',jnfc,'Interpreter','none');
        end
        
        if PLOT_JOINT_V
            jv_uv = uv(joint_i.v);          % joint directional velocity unit vector
            jv_len = jvar*norm(joint_i.v);  % joint directional velocity arrow length
            p1 = p;
            p2 = p + jv_len*jv_uv;
            fv = get_arrow_3d(p1,p2,'color',jvfc,'stemWidth',jvsw,'tipWidth',jvtw,...
                'facealpha',jvfa);
            h{fig_idx,subfig_idx}.jv_arrow{i_idx} = patch(fv,'facecolor',jvfc,'edgeColor','none',...
                'FaceAlpha',jvfa,'FaceLighting','gouraud');
        end
        
        if PLOT_JOINT_W
            jw_uv = uv(joint_i.w);          % joint angular velocity unit vector
            jw_len = jwar*norm(joint_i.w);  % joint angular velocity arrow length
            p1 = p;
            p2 = p + jw_len*jw_uv;
            fv = get_arrow_3d(p1,p2,'color',jwfc,'stemWidth',jwsw,'tipWidth',jwtw,...
                'facealpha',jwfa);
            h{fig_idx,subfig_idx}.jw_arrow{i_idx} = patch(fv,'facecolor',jwfc,'edgeColor','none',...
                'FaceAlpha',jwfa,'FaceLighting','gouraud');
        end
        
    end
    
else
    % -------------------------------------------------------------------------------------------- %
    
    % Plot link (>1)
    if PLOT_LINK
        for i_idx = 1:chain.n_joint
            joint_i = chain.joint(i_idx);
            parent = joint_i.parent;
            if ~isempty(parent)
                joint_fr = chain.joint(parent);
                p1 = joint_fr.p;
                p2 = joint_i.p;
                h{fig_idx,subfig_idx}.line{i_idx}.XData = [p1(1),p2(1)];
                h{fig_idx,subfig_idx}.line{i_idx}.YData = [p1(2),p2(2)];
                h{fig_idx,subfig_idx}.line{i_idx}.ZData = [p1(3),p2(3)];
                
                h{fig_idx,subfig_idx}.line{i_idx}.LineWidth = llw*dragzoom_rate;
            end
        end
    end
    
    % Plot mesh (>1)
    if PLOT_MESH && isfield(chain,'link') % only if link exists
        for i_idx = 1:chain.n_link
            link_i = chain.link(i_idx);
            fv = link_i.fv;
            joint_idx = link_i.joint_idx;
            if (~isempty(joint_idx))
                p_link = chain.joint(joint_idx).p;
                R_link = chain.joint(joint_idx).R;
                % Plot mesh
                if ~isempty(fv) % if mesh exists
                    % Move
                    tform = pr2t(p_link,R_link);
                    set(h{fig_idx,subfig_idx}.mesh_t{i_idx},'Matrix',tform);
                end
            end
        end
    end
    
    if PLOT_BOX && isfield(chain,'link') % (>1)
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            box_i = link_i.box;
            if ~isempty(box_i) && ~isempty(link_i.joint_idx) % if box exists
                joint_idx = link_i.joint_idx;
                p_link = chain.joint(joint_idx).p;
                R_link = chain.joint(joint_idx).R;
                tform = pr2t(p_link,R_link);
                set(h{fig_idx,subfig_idx}.box_t{i_idx},'Matrix',tform);
            end
        end % for i_idx = 1:chain.n_link % for all links
    end
    
    % Plot box added (>1)
    if PLOT_BOX_ADDED && isfield(chain,'link')
        for i_idx = 1:chain.n_link
            link_i = chain.link(i_idx);
            if isfield(link_i,'box_added')
                if ~isempty(link_i.box_added)
                    box_added = link_i.box_added;
                    T = pr2t(chain.joint(link_i.joint_idx).p, chain.joint(link_i.joint_idx).R);
                    if isfield(box_added,'p_offset')
                        p_offset = box_added.p_offset;
                        R_offset = box_added.R_offset;
                        T_offset = pr2t(p_offset,R_offset);
                    else
                        T_offset = pr2t(cv([0,0,0]),eye(3,3));
                    end
                    tform = T*T_offset;
                    set(h{fig_idx,subfig_idx}.box_added_t{i_idx},'Matrix',tform);
                end
            end
        end
    end
    
    % Plot capsule (>1)
    if PLOT_CAPSULE && isfield(chain,'link')
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.capsule)
                cap = link_i.capsule;
                joint_idx = link_i.joint_idx;
                if isempty(joint_idx)
                    joint_idx = get_topmost_idx(chain);
                end
                p_link = chain.joint(joint_idx).p;
                R_link = chain.joint(joint_idx).R;
                
                tform = pr2t(p_link,R_link);
                
                set(h{fig_idx,subfig_idx}.cap_patch_t{i_idx},'Matrix',tform);
            end
        end % for i_idx = 1:chain.n_link % for all links
    end % if PLOT_BOX_ADDED && isfield(chain,'link')
    
    % Plot CoM (>1)
    if PLOT_COM && isfield(chain,'link')
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.com_bar) && ~isempty(link_i.joint_idx)
                joint_idx = link_i.joint_idx;
                com_bar = link_i.com_bar;
                p_joint = chain.joint(joint_idx).p;
                R_joint = chain.joint(joint_idx).R;
                T_joint = pr2t(p_joint,R_joint);
                tform = T_joint*pr2t(com_bar,'');
                set(h{fig_idx,subfig_idx}.com_t{i_idx},'Matrix',tform);
            end
        end % for i_idx = 1:chain.n_link % for all links
    end % if PLOT_COM && isfield(chain,'link')
    
    % Plot link directional velocity (>1)
    if PLOT_LINK_V
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.com_bar)
                joint_idx = link_i.joint_idx;
                com_bar = link_i.com_bar; % in local coordinates (post-multiply this)
                p_joint = chain.joint(joint_idx).p;
                R_joint = chain.joint(joint_idx).R;
                T_joint = pr2t(p_joint,R_joint);
                T_com = T_joint*pr2t(com_bar,'');
                p_com = t2p(T_com);
                lv_uv = uv(link_i.v);          % link directional velocity unit vector
                lv_len = lvar*norm(link_i.v);  % link directional velocity arrow length
                p1 = p_com;
                p2 = p_com + lv_len*lv_uv;
                fv = get_arrow_3d(p1,p2,'color',lvfc,'stemWidth',lvsw,'tipWidth',lvtw,...
                    'facealpha',lvfa);
                h{fig_idx,subfig_idx}.lv_arrow{i_idx}.Faces     = fv.faces;
                h{fig_idx,subfig_idx}.lv_arrow{i_idx}.Vertices  = fv.vertices;
            end
        end
    end
    
    % Plot link angular velocity (>1)
    if PLOT_LINK_W
        for i_idx = 1:chain.n_link % for all links
            link_i = chain.link(i_idx);
            if ~isempty(link_i.com_bar)
                joint_idx = link_i.joint_idx;
                com_bar = link_i.com_bar; % in local coordinates (post-multiply this)
                p_joint = chain.joint(joint_idx).p;
                R_joint = chain.joint(joint_idx).R;
                T_joint = pr2t(p_joint,R_joint);
                T_com = T_joint*pr2t(com_bar,'');
                p_com = t2p(T_com);
                lw_uv = uv(link_i.w);          % link directional velocity unit vector
                lw_len = lwar*norm(link_i.w);  % link directional velocity arrow length
                p1 = p_com;
                p2 = p_com + lw_len*lw_uv;
                fv = get_arrow_3d(p1,p2,'color',lwfc,'stemWidth',lwsw,'tipWidth',lwtw,...
                    'facealpha',lwfa);
                h{fig_idx,subfig_idx}.lw_arrow{i_idx}.Faces     = fv.faces;
                h{fig_idx,subfig_idx}.lw_arrow{i_idx}.Vertices  = fv.vertices;
            end
        end
    end
    
    % Plot joint related information (>1)
    for i_idx = 1:chain.n_joint % plot joint information
        joint_i = chain.joint(i_idx);
        p = joint_i.p;
        R = joint_i.R;
        ex = R(:,1); ey = R(:,2); ez = R(:,3);
        
        if PLOT_JOINT_AXIS
            h{fig_idx,subfig_idx}.cx{i_idx}.XData = [p(1),p(1)+jal*ex(1)];
            h{fig_idx,subfig_idx}.cx{i_idx}.YData = [p(2),p(2)+jal*ex(2)];
            h{fig_idx,subfig_idx}.cx{i_idx}.ZData = [p(3),p(3)+jal*ex(3)];
            h{fig_idx,subfig_idx}.cy{i_idx}.XData = [p(1),p(1)+jal*ey(1)];
            h{fig_idx,subfig_idx}.cy{i_idx}.YData = [p(2),p(2)+jal*ey(2)];
            h{fig_idx,subfig_idx}.cy{i_idx}.ZData = [p(3),p(3)+jal*ey(3)];
            h{fig_idx,subfig_idx}.cz{i_idx}.XData = [p(1),p(1)+jal*ez(1)];
            h{fig_idx,subfig_idx}.cz{i_idx}.YData = [p(2),p(2)+jal*ez(2)];
            h{fig_idx,subfig_idx}.cz{i_idx}.ZData = [p(3),p(3)+jal*ez(3)];
            h{fig_idx,subfig_idx}.cx{i_idx}.LineWidth = jalw*dragzoom_rate;
            h{fig_idx,subfig_idx}.cy{i_idx}.LineWidth = jalw*dragzoom_rate;
            h{fig_idx,subfig_idx}.cz{i_idx}.LineWidth = jalw*dragzoom_rate;
        end
        
        if PLOT_JOINT_SPHERE
            tform = pr2t(p,'');
            set(h{fig_idx,subfig_idx}.sphere_t{i_idx},'Matrix',tform);
        end
        
        if PLOT_ROTATE_AXIS
            a = joint_i.a; % axis
            if sum(abs(a)) ~= 0
                p1 = p;
                p2 = p+R*a*ral;
                fv = get_arrow_3d(p1,p2,'color','b','stemWidth',rasw,'tipWidth',ratw,...
                    'facealpha',raa);
                h{fig_idx,subfig_idx}.arrow{i_idx}.Faces = fv.faces;
                h{fig_idx,subfig_idx}.arrow{i_idx}.Vertices = fv.vertices;
            end
        end
        
        if PLOT_JOINT_NAME || PLOT_JOINT_TORQUE
            text_p = p;
            if PLOT_JOINT_TORQUE && isequal(joint_i.type,'revolute')
                text_str = sprintf(' [%d]%s(%.2f)',i_idx,joint_i.name,joint_i.u);
            else
                text_str = sprintf(' [%d]%s',i_idx,joint_i.name);
            end
            h{fig_idx,subfig_idx}.text{i_idx}.Position = text_p;
            h{fig_idx,subfig_idx}.text{i_idx}.String   = text_str;
            h{fig_idx,subfig_idx}.text{i_idx}.FontSize = dragzoom_rate*jnfs; % adaptive font size
        end
        
        if PLOT_JOINT_V
            
            jv_uv = uv(joint_i.v);          % joint directional velocity unit vector
            jv_len = jvar*norm(joint_i.v);  % joint directional velocity arrow length
            p1 = p;
            p2 = p + jv_len*jv_uv;
            fv = get_arrow_3d(p1,p2,'color',jvfc,'stemWidth',jvsw,'tipWidth',jvtw,...
                'facealpha',jvfa);
            h{fig_idx,subfig_idx}.jv_arrow{i_idx}.Faces = fv.faces;
            h{fig_idx,subfig_idx}.jv_arrow{i_idx}.Vertices = fv.vertices;
            
        end
        
        if PLOT_JOINT_W
            
            jw_uv = uv(joint_i.w);          % joint angular velocity unit vector
            jw_len = jwar*norm(joint_i.w);  % joint angular velocity arrow length
            p1 = p;
            p2 = p + jw_len*jw_uv;
            fv = get_arrow_3d(p1,p2,'color',jwfc,'stemWidth',jwsw,'tipWidth',jwtw,...
                'facealpha',jwfa);
            h{fig_idx,subfig_idx}.jw_arrow{i_idx}.Faces = fv.faces;
            h{fig_idx,subfig_idx}.jw_arrow{i_idx}.Vertices = fv.vertices;
            
        end
    end % for i_idx = 1:chain.n_joint % plot joint information
    
end

fig = h{fig_idx,subfig_idx}.fig;
