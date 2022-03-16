function [fig,h_coord,h_sphere] = plot_T(T,varargin)
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
addParameter(ps,'PLOT_AXIS',true);       % plot axis 
addParameter(ps,'all',1.0);              % axis line length
addParameter(ps,'alc','');               % axis line color
addParameter(ps,'alw',2);                % axis line width
addParameter(ps,'als','-');              % axis line style
addParameter(ps,'PLOT_AXIS_TIP',false);  % plot axis tip
addParameter(ps,'atr',0.05);             % axis tip rate w.r.t. axis line length
addParameter(ps,'PLOT_SPHERE',false);    % plot sphere at the center
addParameter(ps,'sr',0.1);               % sphere radius
addParameter(ps,'sfc',0.5*[1,1,1]);      % sphere face color
addParameter(ps,'sfa',0.5);              % sphere face alpha
addParameter(ps,'text_str','');          % show text 
addParameter(ps,'text_fs',15);           % text font size
addParameter(ps,'text_fn','consolas');   % text font name
addParameter(ps,'text_color','k');       % text color 
addParameter(ps,'text_interp','none');   % text interpreter
addParameter(ps,'TEXT_AT_ZTIP',0);
addParameter(ps,'USE_ZOOMRATE',1);
parse(ps,varargin{:});
fig_idx         = ps.Results.fig_idx;
subfig_idx      = ps.Results.subfig_idx;
PLOT_AXIS       = ps.Results.PLOT_AXIS;
all             = ps.Results.all;
alc             = ps.Results.alc;
alw             = ps.Results.alw;
als             = ps.Results.als;
PLOT_AXIS_TIP   = ps.Results.PLOT_AXIS_TIP;
atr             = ps.Results.atr;
PLOT_SPHERE     = ps.Results.PLOT_SPHERE;
sr              = ps.Results.sr;
sfc             = ps.Results.sfc;
sfa             = ps.Results.sfa;
text_str        = ps.Results.text_str;
text_fs         = ps.Results.text_fs;
text_fn         = ps.Results.text_fn;
text_color      = ps.Results.text_color;
text_interp     = ps.Results.text_interp;
TEXT_AT_ZTIP    = ps.Results.TEXT_AT_ZTIP;
USE_ZOOMRATE    = ps.Results.USE_ZOOMRATE;

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

% Get p and R of T
p = T([1,2,3],4);
R = T([1,2,3],[1,2,3]);
ex = R(:,1); ey = R(:,2); ez = R(:,3);

% Plot start here
if h{fig_idx,subfig_idx}.first_flag || ~ishandle(h{fig_idx,subfig_idx}.fig)
    h{fig_idx,subfig_idx}.first_flag = false;
    h{fig_idx,subfig_idx}.fig = figure(fig_idx);
    hold on;
    if PLOT_AXIS
        if isempty(alc)
            h{fig_idx,subfig_idx}.cx = plot3([p(1),p(1)+all*ex(1)],[p(2),p(2)+all*ex(2)],...
                [p(3),p(3)+all*ex(3)],'r','LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.cy = plot3([p(1),p(1)+all*ey(1)],[p(2),p(2)+all*ey(2)],...
                [p(3),p(3)+all*ey(3)],'g','LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.cz = plot3([p(1),p(1)+all*ez(1)],[p(2),p(2)+all*ez(2)],...
                [p(3),p(3)+all*ez(3)],'b','LineWidth',alw*dragzoom_rate,'LineStyle',als);
        else
            h{fig_idx,subfig_idx}.cx = plot3([p(1),p(1)+all*ex(1)],[p(2),p(2)+all*ex(2)],...
                [p(3),p(3)+all*ex(3)],alc,'LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.cy = plot3([p(1),p(1)+all*ey(1)],[p(2),p(2)+all*ey(2)],...
                [p(3),p(3)+all*ey(3)],alc,'LineWidth',alw*dragzoom_rate,'LineStyle',als);
            h{fig_idx,subfig_idx}.cz = plot3([p(1),p(1)+all*ez(1)],[p(2),p(2)+all*ez(2)],...
                [p(3),p(3)+all*ez(3)],alc,'LineWidth',alw*dragzoom_rate,'LineStyle',als);
        end
    end
    
    if PLOT_AXIS_TIP
        [xs,ys,zs] = sphere(20); 
        fv = surf2patch(atr*all*xs,atr*all*ys,atr*all*zs);
        
        ax = p(1)+all*ex(1);
        ay = p(2)+all*ex(2);
        az = p(3)+all*ex(3);
        h{fig_idx,subfig_idx}.sx = patch(fv,...
            'EdgeColor','none','FaceColor','r','FaceAlpha',0.9,'facelighting','gouraud');
        h{fig_idx,subfig_idx}.sx_t = hgtransform;
        set(h{fig_idx,subfig_idx}.sx,'parent',h{fig_idx,subfig_idx}.sx_t);
        tform = pr2t([ax,ay,az],eye(3,3));
        set(h{fig_idx,subfig_idx}.sx_t,'Matrix',tform);
        
        ax = p(1)+all*ey(1);
        ay = p(2)+all*ey(2);
        az = p(3)+all*ey(3);
        h{fig_idx,subfig_idx}.sy = patch(fv,...
            'EdgeColor','none','FaceColor','g','FaceAlpha',0.9,'facelighting','gouraud');
        h{fig_idx,subfig_idx}.sy_t = hgtransform;
        set(h{fig_idx,subfig_idx}.sy,'parent',h{fig_idx,subfig_idx}.sy_t);
        tform = pr2t([ax,ay,az],eye(3,3));
        set(h{fig_idx,subfig_idx}.sy_t,'Matrix',tform);
        
        ax = p(1)+all*ez(1);
        ay = p(2)+all*ez(2);
        az = p(3)+all*ez(3);
        h{fig_idx,subfig_idx}.sz = patch(fv,...
            'EdgeColor','none','FaceColor','b','FaceAlpha',0.9,'facelighting','gouraud');
        h{fig_idx,subfig_idx}.sz_t = hgtransform;
        set(h{fig_idx,subfig_idx}.sz,'parent',h{fig_idx,subfig_idx}.sz_t);
        tform = pr2t([ax,ay,az],eye(3,3));
        set(h{fig_idx,subfig_idx}.sz_t,'Matrix',tform);
    end
    
    if PLOT_SPHERE
        [x,y,z] = ellipsoid(0,0,0,sr,sr,sr,30);
        fv = surf2patch(x,y,z);
        h{fig_idx,subfig_idx}.sphere = patch(fv,...
            'EdgeColor','none','FaceColor',sfc,'FaceAlpha',sfa,'facelighting','gouraud');
        h{fig_idx,subfig_idx}.sphere_t = hgtransform;
        set(h{fig_idx,subfig_idx}.sphere,'parent',h{fig_idx,subfig_idx}.sphere_t);
        tform = pr2t(p,eye(3,3));
        set(h{fig_idx,subfig_idx}.sphere_t,'Matrix',tform);
    end
    
    if ~isempty(text_str)
        text_p = p;
        if TEXT_AT_ZTIP
            text_p = text_p + all*ez;
        end
        h{fig_idx,subfig_idx}.text = text(text_p(1),text_p(2),text_p(3),[' ',text_str],...
            'FontSize',text_fs*dragzoom_rate,'FontName',text_fn,'Color',text_color,...
            'Interpreter',text_interp);
    end
    
else
    
    if PLOT_AXIS
        h{fig_idx,subfig_idx}.cx.XData = [p(1),p(1)+all*ex(1)];
        h{fig_idx,subfig_idx}.cx.YData = [p(2),p(2)+all*ex(2)];
        h{fig_idx,subfig_idx}.cx.ZData = [p(3),p(3)+all*ex(3)];
        h{fig_idx,subfig_idx}.cy.XData = [p(1),p(1)+all*ey(1)];
        h{fig_idx,subfig_idx}.cy.YData = [p(2),p(2)+all*ey(2)];
        h{fig_idx,subfig_idx}.cy.ZData = [p(3),p(3)+all*ey(3)];
        h{fig_idx,subfig_idx}.cz.XData = [p(1),p(1)+all*ez(1)];
        h{fig_idx,subfig_idx}.cz.YData = [p(2),p(2)+all*ez(2)];
        h{fig_idx,subfig_idx}.cz.ZData = [p(3),p(3)+all*ez(3)];
        
        h{fig_idx,subfig_idx}.cx.LineWidth = dragzoom_rate*alw;
        h{fig_idx,subfig_idx}.cy.LineWidth = dragzoom_rate*alw;
        h{fig_idx,subfig_idx}.cz.LineWidth = dragzoom_rate*alw;
    end
    
    if PLOT_AXIS_TIP
        ax = p(1)+all*ex(1);
        ay = p(2)+all*ex(2);
        az = p(3)+all*ex(3);
        tform = pr2t([ax,ay,az],eye(3,3));
        set(h{fig_idx,subfig_idx}.sx_t,'Matrix',tform);
        
        ax = p(1)+all*ey(1);
        ay = p(2)+all*ey(2);
        az = p(3)+all*ey(3);
        tform = pr2t([ax,ay,az],eye(3,3));
        set(h{fig_idx,subfig_idx}.sy_t,'Matrix',tform);
        
        ax = p(1)+all*ez(1);
        ay = p(2)+all*ez(2);
        az = p(3)+all*ez(3);
        tform = pr2t([ax,ay,az],eye(3,3));
        set(h{fig_idx,subfig_idx}.sz_t,'Matrix',tform);
    end
    
    if PLOT_SPHERE
        tform = pr2t(p,eye(3,3));
        set(h{fig_idx,subfig_idx}.sphere_t,'Matrix',tform);
    end
    
    if ~isempty(text_str)
        text_p = p;
        if TEXT_AT_ZTIP
            text_p = text_p + all*ez;
        end
        h{fig_idx,subfig_idx}.text.Position = text_p;
        h{fig_idx,subfig_idx}.text.FontSize = dragzoom_rate*text_fs;
        h{fig_idx,subfig_idx}.text.String = [' ',text_str];
    end
    
end

% Output arguments 
fig = h{fig_idx,subfig_idx}.fig;
if PLOT_AXIS
    h_coord = h{fig_idx,subfig_idx}.cx;
else
    h_coord = '';
end
if PLOT_SPHERE
    h_sphere= h{fig_idx,subfig_idx}.sphere;
else
    h_sphere = '';   
end
