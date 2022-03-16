function plot_line(p1,p2,varargin)
%
% Plot a straight line between two points p1 and p2 in 3D
%
global zoomPct % zoom rate from dragzoom 
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:10, for j = 1:100, h{i,j}.first_flag = true; end; end; end

% Parse options
ps = inputParser;
addParameter(ps,'fig_idx',1);
addParameter(ps,'subfig_idx',1);
addParameter(ps,'lc','k');               % line color 
addParameter(ps,'lw',2);                 % line width
addParameter(ps,'ls','-');               % line style
addParameter(ps,'text_str','');
addParameter(ps,'text_fs',13);
addParameter(ps,'text_fn','consolas');
addParameter(ps,'text_color','k');
addParameter(ps,'PLOT_LINE_TIP',0);
addParameter(ps,'ltfc','r');             % line tip face color
addParameter(ps,'ltfa',0.5);             % line tip face alpha
addParameter(ps,'ltr',0.05);             % line tip radius
addParameter(ps,'USE_ZOOMRATE',1);
parse(ps,varargin{:});
fig_idx         = ps.Results.fig_idx;
subfig_idx      = ps.Results.subfig_idx;
lc              = ps.Results.lc;
lw              = ps.Results.lw;
ls              = ps.Results.ls;
text_str        = ps.Results.text_str;
text_fs         = ps.Results.text_fs;
text_fn         = ps.Results.text_fn;
text_color      = ps.Results.text_color;
PLOT_LINE_TIP   = ps.Results.PLOT_LINE_TIP;
ltfc            = ps.Results.ltfc;
ltfa            = ps.Results.ltfa;
ltr             = ps.Results.ltr;
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

% Plot start here
if h{fig_idx,subfig_idx}.first_flag || ~ishandle(h{fig_idx,subfig_idx}.fig)
    h{fig_idx,subfig_idx}.first_flag = false;
    h{fig_idx,subfig_idx}.fig = figure(fig_idx);
    hold on;
    
    h{fig_idx,subfig_idx}.line = plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],...
        'Color',lc,'LineWidth',lw*dragzoom_rate,'LineStyle',ls);
    
    if ~isempty(text_str)
        p = 0.5*(p1+p2);
        h{fig_idx,subfig_idx}.text = text(p(1),p(2),p(3),[' ',text_str],...
            'FontSize',text_fs*dragzoom_rate,'FontName',text_fn,'Color',text_color,...
            'Interpreter','none');
    end
    
    if PLOT_LINE_TIP
        [xs,ys,zs] = sphere(20); 
        fv = surf2patch(ltr*xs,ltr*ys,ltr*zs);
        
        h{fig_idx,subfig_idx}.s1 = patch(fv,...
            'EdgeColor','none','FaceColor',ltfc,'FaceAlpha',ltfa,'facelighting','gouraud');
        h{fig_idx,subfig_idx}.s1_t = hgtransform;
        set(h{fig_idx,subfig_idx}.s1,'parent',h{fig_idx,subfig_idx}.s1_t);
        tform = pr2t(cv(p1),eye(3,3));
        set(h{fig_idx,subfig_idx}.s1_t,'Matrix',tform);
        
        h{fig_idx,subfig_idx}.s2 = patch(fv,...
            'EdgeColor','none','FaceColor',ltfc,'FaceAlpha',ltfa,'facelighting','gouraud');
        h{fig_idx,subfig_idx}.s2_t = hgtransform;
        set(h{fig_idx,subfig_idx}.s2,'parent',h{fig_idx,subfig_idx}.s2_t);
        tform = pr2t(cv(p2),eye(3,3));
        set(h{fig_idx,subfig_idx}.s2_t,'Matrix',tform);
    end
    
else
    
    h{fig_idx,subfig_idx}.line.XData = [p1(1),p2(1)];
    h{fig_idx,subfig_idx}.line.YData = [p1(2),p2(2)];
    h{fig_idx,subfig_idx}.line.ZData = [p1(3),p2(3)];
    h{fig_idx,subfig_idx}.line.LineWidth = lw*dragzoom_rate;
    
    if ~isempty(text_str)
        p = 0.5*(p1+p2);
        h{fig_idx,subfig_idx}.text.Position = p;
        h{fig_idx,subfig_idx}.text.FontSize = text_fs*dragzoom_rate;
        h{fig_idx,subfig_idx}.text.String = [' ',text_str];
    end
    
    if PLOT_LINE_TIP
        tform = pr2t(cv(p1),eye(3,3));
        set(h{fig_idx,subfig_idx}.s1_t,'Matrix',tform);
        
        tform = pr2t(cv(p2),eye(3,3));
        set(h{fig_idx,subfig_idx}.s2_t,'Matrix',tform);
    end
    
end
