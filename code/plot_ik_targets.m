function plot_ik_targets(varargin)
%
% Plot homogeneous transformation Matrix
%
global zoomPct % zoom rate from dragzoom
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:10, for j = 1:100, h{i,j}.first_flag = true; end; end; end

% Parse options
ps = inputParser;
addParameter(ps,'fig_idx',1);
addParameter(ps,'subfig_idx',1);
addParameter(ps,'chain_robot','');              % robot chain
addParameter(ps,'ik_plot_info','');             % ik plot information
addParameter(ps,'PLOT_AXIS',true);              % plot axis
addParameter(ps,'all',1.0);                     % axis line length
addParameter(ps,'alw',2);                       % axis line width
addParameter(ps,'als','-');                     % axis line style
addParameter(ps,'PLOT_SPHERE',true);            % plot sphere at the center
addParameter(ps,'sr',0.1);                      % sphere radius
addParameter(ps,'sfc','');                      % sphere face color
addParameter(ps,'sfa',0.5);                     % sphere face alpha
addParameter(ps,'PLOT_ARROW',false);            % plot direction arrow
addParameter(ps,'adl',0.2);                     % arrow direction length
addParameter(ps,'adc','r');                     % arrow direction color
addParameter(ps,'adsw',0.02);                   % arrow direction stem width
addParameter(ps,'adtw',0.04);                   % arrow direction tip width
addParameter(ps,'USE_ZOOMRATE',0);
addParameter(ps,'RESET',0);
parse(ps,varargin{:});
fig_idx                 = ps.Results.fig_idx;
subfig_idx              = ps.Results.subfig_idx;
chain_robot             = ps.Results.chain_robot;
ik_plot_info            = ps.Results.ik_plot_info;
PLOT_AXIS               = ps.Results.PLOT_AXIS;
all                     = ps.Results.all;
alw                     = ps.Results.alw;
als                     = ps.Results.als;
PLOT_SPHERE             = ps.Results.PLOT_SPHERE;
sr                      = ps.Results.sr;
sfc                     = ps.Results.sfc;
sfa                     = ps.Results.sfa;
PLOT_ARROW              = ps.Results.PLOT_ARROW;
adl                     = ps.Results.adl;
adc                     = ps.Results.adc;
adsw                    = ps.Results.adsw;
adtw                    = ps.Results.adtw;
USE_ZOOMRATE            = ps.Results.USE_ZOOMRATE;
RESET                   = ps.Results.RESET;

% Dragzoom rate
if isempty(zoomPct) || (USE_ZOOMRATE==0)
    dragzoom_rate = 1;
else
    dragzoom_rate = 1 + (zoomPct-100)/150; % normalize by 100
    max_dragzoom_rate = 5;
    if dragzoom_rate > max_dragzoom_rate
        dragzoom_rate = max_dragzoom_rate;
    elseif dragzoom_rate < 1/max_dragzoom_rate
        dragzoom_rate = 1/max_dragzoom_rate;
    end
end

% Reset arrow direction
if PLOT_ARROW || RESET
    if isfield(h{fig_idx,subfig_idx},'arrow')
        for i_idx = 1:length(h{fig_idx,subfig_idx}.arrow)
            delete(h{fig_idx,subfig_idx}.arrow{i_idx});
        end
    end
end

% Reset all the handles 
if RESET
    for fig_idx = 1:10
        for subfig_idx = 1:100
            if ~h{fig_idx,subfig_idx}.first_flag
                
                if isfield(h{fig_idx,subfig_idx},'ik_p_curr')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ik_p_curr)
                        delete(h{fig_idx,subfig_idx}.ik_p_curr{i_idx});
                    end
                end
                
                if isfield(h{fig_idx,subfig_idx},'ik_p_trgt')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ik_p_trgt)
                        delete(h{fig_idx,subfig_idx}.ik_p_trgt{i_idx});
                    end
                end
                
                if isfield(h{fig_idx,subfig_idx},'ex_curr')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ex_curr)
                        delete(h{fig_idx,subfig_idx}.ex_curr{i_idx});
                    end
                end
                
                if isfield(h{fig_idx,subfig_idx},'ey_curr')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ey_curr)
                        delete(h{fig_idx,subfig_idx}.ey_curr{i_idx});
                    end
                end
                
                if isfield(h{fig_idx,subfig_idx},'ez_curr')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ez_curr)
                        delete(h{fig_idx,subfig_idx}.ez_curr{i_idx});
                    end
                end
                
                if isfield(h{fig_idx,subfig_idx},'ex_trgt')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ex_trgt)
                        delete(h{fig_idx,subfig_idx}.ex_trgt{i_idx});
                    end
                end
                
                if isfield(h{fig_idx,subfig_idx},'ey_trgt')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ey_trgt)
                        delete(h{fig_idx,subfig_idx}.ey_trgt{i_idx});
                    end
                end
                
                if isfield(h{fig_idx,subfig_idx},'ez_trgt')
                    for i_idx = 1:length(h{fig_idx,subfig_idx}.ez_trgt)
                        delete(h{fig_idx,subfig_idx}.ez_trgt{i_idx});
                    end
                end
                
                h{fig_idx,subfig_idx}.first_flag = true; % revert flag
            end
        end
    end
    
    return;
end

% Following information are needed to plot 
ik_joint_names      = ik_plot_info.joint_names;
ik_types            = ik_plot_info.types;
ik_trgt_coords      = ik_plot_info.trgt_coords;

% Plot start here
if h{fig_idx,subfig_idx}.first_flag || ~ishandle(h{fig_idx,subfig_idx}.fig)
    h{fig_idx,subfig_idx}.first_flag = false;
    h{fig_idx,subfig_idx}.fig = figure(fig_idx);
    colors = linspecer(length(ik_joint_names));
    [x,y,z] = ellipsoid(0,0,0,sr,sr,sr,30);
    fv_sphere = surf2patch(x,y,z);
    for i_idx = 1:length(ik_joint_names) % for different IK targets
        ik_joint_name = ik_joint_names{i_idx};  % ik joint name
        ik_type       = ik_types{i_idx};        % ik target type
        ik_trgt_coord = ik_trgt_coords{i_idx};  % ik target coordinates
        
        IK_P = 0; IK_R = 0;
        switch ik_type
            case 'IK_P'
                IK_P = 1;
            case 'IK_R'
                IK_R = 1;
            case 'IK_PR'
                IK_P = 1;
                IK_R = 1;
            otherwise
                fprintf(2,'[plot_ik_targets] Unknown ik_type:[%s].\n',ik_type);
        end
        % Current position and rotation of the robot
        joint_idx = idx_cell(chain_robot.joint_names,ik_joint_name);
        p_curr = chain_robot.joint(joint_idx).p;
        R_curr = chain_robot.joint(joint_idx).R;
        
        % IK target position and rotation
        [p_trgt,R_trgt] = t2pr(ik_trgt_coord);
        
        % Plot sphere (for position targets)
        if PLOT_SPHERE && IK_P
            if isempty(sfc)
                sfc_i = colors(i_idx,:);
            else
                sfc_i = sfc;
            end
            h{fig_idx,subfig_idx}.ik_p_curr{i_idx} = patch(fv_sphere,...
                'EdgeColor','none','FaceColor',sfc_i,'FaceAlpha',sfa,'facelighting','gouraud');
            h{fig_idx,subfig_idx}.ik_p_curr_t{i_idx} = hgtransform;
            set(h{fig_idx,subfig_idx}.ik_p_curr{i_idx},...
                'parent',h{fig_idx,subfig_idx}.ik_p_curr_t{i_idx});
            tform = pr2t(p_curr,eye(3,3));
            set(h{fig_idx,subfig_idx}.ik_p_curr_t{i_idx},'Matrix',tform);
            
            h{fig_idx,subfig_idx}.ik_p_trgt{i_idx} = patch(fv_sphere,...
                'EdgeColor','none','FaceColor',sfc_i,'FaceAlpha',sfa,'facelighting','gouraud');
            h{fig_idx,subfig_idx}.ik_p_trgt_t{i_idx} = hgtransform;
            set(h{fig_idx,subfig_idx}.ik_p_trgt{i_idx},...
                'parent',h{fig_idx,subfig_idx}.ik_p_trgt_t{i_idx});
            tform = pr2t(p_trgt,eye(3,3));
            set(h{fig_idx,subfig_idx}.ik_p_trgt_t{i_idx},'Matrix',tform);
        end
        
        % Plot arrow direction 
        if PLOT_ARROW && IK_P
            p_curr2trgt = p_trgt - p_curr;
            p1 = p_curr;
            p2 = p1 + adl*uv(p_curr2trgt);
            fv_arrow = get_arrow_3d(p1,p2,...
                'color',adc,'stemWidth',adsw,'tipWidth',adtw,'facealpha',0.5);
            h{fig_idx,subfig_idx}.arrow{i_idx} = patch(fv_arrow,...
                'facecolor',adc,'edgeColor','none','FaceAlpha',0.5,'FaceLighting','gouraud');
        end
        
        % Plot axis (for rotation targets)
        if PLOT_AXIS && IK_R
            p = p_curr;
            R = R_curr;
            ex = R(:,1); ey = R(:,2); ez = R(:,3);
            h{fig_idx,subfig_idx}.ex_curr{i_idx} = plot3([p(1),p(1)+all*ex(1)],[p(2),p(2)+all*ex(2)],...
                [p(3),p(3)+all*ex(3)],'r','LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.ey_curr{i_idx} = plot3([p(1),p(1)+all*ey(1)],[p(2),p(2)+all*ey(2)],...
                [p(3),p(3)+all*ey(3)],'g','LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.ez_curr{i_idx} = plot3([p(1),p(1)+all*ez(1)],[p(2),p(2)+all*ez(2)],...
                [p(3),p(3)+all*ez(3)],'b','LineWidth',alw*dragzoom_rate,'LineStyle',als);
            
            p = p_curr;
            R = R_trgt;
            ex = R(:,1); ey = R(:,2); ez = R(:,3);
            h{fig_idx,subfig_idx}.ex_trgt{i_idx} = plot3([p(1),p(1)+all*ex(1)],[p(2),p(2)+all*ex(2)],...
                [p(3),p(3)+all*ex(3)],'r','LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.ey_trgt{i_idx} = plot3([p(1),p(1)+all*ey(1)],[p(2),p(2)+all*ey(2)],...
                [p(3),p(3)+all*ey(3)],'g','LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.ez_trgt{i_idx} = plot3([p(1),p(1)+all*ez(1)],[p(2),p(2)+all*ez(2)],...
                [p(3),p(3)+all*ez(3)],'b','LineWidth',alw*dragzoom_rate,'LineStyle',als);
        end
        
    end % for i_idx = 1:length(ik_joint_names) % for different IK targets
    
else
    
    for i_idx = 1:length(ik_joint_names) % for different IK targets
        ik_joint_name = ik_joint_names{i_idx};  % ik joint name
        ik_type       = ik_types{i_idx};        % ik target type
        ik_trgt_coord = ik_trgt_coords{i_idx};  % ik target coordinates
        IK_P = 0; IK_R = 0;
        switch ik_type
            case 'IK_P'
                IK_P = 1;
            case 'IK_R'
                IK_R = 1;
            case 'IK_PR'
                IK_P = 1;
                IK_R = 1;
            otherwise
                fprintf(2,'[plot_ik_targets] Unknown ik_type:[%s].\n',ik_type);
        end
        % Current position and rotation
        joint_idx = idx_cell(chain_robot.joint_names,ik_joint_name);
        p_curr = chain_robot.joint(joint_idx).p;
        R_curr = chain_robot.joint(joint_idx).R;
        
        % IK target position and rotation
        [p_trgt,R_trgt] = t2pr(ik_trgt_coord);
        
        % Plot sphere (for position targets)
        if PLOT_SPHERE && IK_P
            tform = pr2t(p_curr,eye(3,3));
            set(h{fig_idx,subfig_idx}.ik_p_curr_t{i_idx},'Matrix',tform);
            
            tform = pr2t(p_trgt,eye(3,3));
            set(h{fig_idx,subfig_idx}.ik_p_trgt_t{i_idx},'Matrix',tform);
        end
        
        % Plot arrow direction 
        if PLOT_ARROW && IK_P
            p_curr2trgt = p_trgt - p_curr;
            p1 = p_curr;
            p2 = p1 + adl*uv(p_curr2trgt);
            fv_arrow = get_arrow_3d(p1,p2,...
                'color',adc,'stemWidth',adsw,'tipWidth',adtw,'facealpha',0.5);
            h{fig_idx,subfig_idx}.arrow{i_idx} = patch(fv_arrow,...
                'facecolor',adc,'edgeColor','none','FaceAlpha',0.5,'FaceLighting','gouraud');
        end
        
        % Plot axis (for rotation targets)
        if PLOT_AXIS && IK_R
            p = p_curr;
            R = R_curr;
            ex = R(:,1); ey = R(:,2); ez = R(:,3);
            h{fig_idx,subfig_idx}.ex_curr{i_idx}.XData = [p(1),p(1)+all*ex(1)];
            h{fig_idx,subfig_idx}.ex_curr{i_idx}.YData = [p(2),p(2)+all*ex(2)];
            h{fig_idx,subfig_idx}.ex_curr{i_idx}.ZData = [p(3),p(3)+all*ex(3)];
            h{fig_idx,subfig_idx}.ex_curr{i_idx}.LineWidth = alw*dragzoom_rate;
            
            h{fig_idx,subfig_idx}.ey_curr{i_idx}.XData = [p(1),p(1)+all*ey(1)];
            h{fig_idx,subfig_idx}.ey_curr{i_idx}.YData = [p(2),p(2)+all*ey(2)];
            h{fig_idx,subfig_idx}.ey_curr{i_idx}.ZData = [p(3),p(3)+all*ey(3)];
            h{fig_idx,subfig_idx}.ey_curr{i_idx}.LineWidth = alw*dragzoom_rate;
            
            h{fig_idx,subfig_idx}.ez_curr{i_idx}.XData = [p(1),p(1)+all*ez(1)];
            h{fig_idx,subfig_idx}.ez_curr{i_idx}.YData = [p(2),p(2)+all*ez(2)];
            h{fig_idx,subfig_idx}.ez_curr{i_idx}.ZData = [p(3),p(3)+all*ez(3)];
            h{fig_idx,subfig_idx}.ez_curr{i_idx}.LineWidth = alw*dragzoom_rate;
            
            p = p_curr;
            R = R_trgt;
            ex = R(:,1); ey = R(:,2); ez = R(:,3);
            h{fig_idx,subfig_idx}.ex_trgt{i_idx}.XData = [p(1),p(1)+all*ex(1)];
            h{fig_idx,subfig_idx}.ex_trgt{i_idx}.YData = [p(2),p(2)+all*ex(2)];
            h{fig_idx,subfig_idx}.ex_trgt{i_idx}.ZData = [p(3),p(3)+all*ex(3)];
            h{fig_idx,subfig_idx}.ex_trgt{i_idx}.LineWidth = alw*dragzoom_rate;
            
            h{fig_idx,subfig_idx}.ey_trgt{i_idx}.XData = [p(1),p(1)+all*ey(1)];
            h{fig_idx,subfig_idx}.ey_trgt{i_idx}.YData = [p(2),p(2)+all*ey(2)];
            h{fig_idx,subfig_idx}.ey_trgt{i_idx}.ZData = [p(3),p(3)+all*ey(3)];
            h{fig_idx,subfig_idx}.ey_trgt{i_idx}.LineWidth = alw*dragzoom_rate;
            
            h{fig_idx,subfig_idx}.ez_trgt{i_idx}.XData = [p(1),p(1)+all*ez(1)];
            h{fig_idx,subfig_idx}.ez_trgt{i_idx}.YData = [p(2),p(2)+all*ez(2)];
            h{fig_idx,subfig_idx}.ez_trgt{i_idx}.ZData = [p(3),p(3)+all*ez(3)];
            h{fig_idx,subfig_idx}.ez_trgt{i_idx}.LineWidth = alw*dragzoom_rate;
        end
        
    end % for i_idx = 1:length(ik_joint_names) % for different IK targets
    
end



