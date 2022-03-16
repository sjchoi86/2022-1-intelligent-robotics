function plot_arrow_3d(p1,p2,varargin)
%
% Plot 3D arrow
%
global zoomPct % zoom rate from dragzoom 
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:10, for j = 1:100, h{i,j}.first_flag = true; end; end; end

% Parse input arguments
if isempty(p1), p1 = cv([0,0,0]); end
len = norm(p1-p2);
ps = inputParser;
addParameter(ps,'fig_idx',1);
addParameter(ps,'subfig_idx',1);
addParameter(ps,'alpha',0.5);
addParameter(ps,'color','r');
addParameter(ps,'sw',len/10);    % stem width
addParameter(ps,'tw',len/5);     % tip width
addParameter(ps,'text_str','');
addParameter(ps,'text_fs',13);
addParameter(ps,'text_fn','consolas');
addParameter(ps,'text_color','k');
addParameter(ps,'interpreter','none');
addParameter(ps,'USE_ZOOMRATE',1);
parse(ps,varargin{:});
fig_idx         = ps.Results.fig_idx;
subfig_idx      = ps.Results.subfig_idx;
alpha           = ps.Results.alpha;
color           = ps.Results.color;
sw              = ps.Results.sw;
tw              = ps.Results.tw;
text_str        = ps.Results.text_str;
text_fs         = ps.Results.text_fs;
text_fn         = ps.Results.text_fn;
text_color      = ps.Results.text_color;
interpreter     = ps.Results.interpreter;
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

% Get 3D arrow mesh
fv = get_arrow_3d(p1,p2,'color',color,'stemWidth',sw,'tipWidth',tw,'facealpha',alpha);
if h{fig_idx,subfig_idx}.first_flag || (~ishandle(h{fig_idx,subfig_idx}.fig))
    h{fig_idx,subfig_idx}.first_flag = false;
    h{fig_idx,subfig_idx}.fig = figure(fig_idx);
    % Plot arrow
    h{fig_idx,subfig_idx}.arrow = patch(fv,'facecolor',color,'edgeColor','none',...
        'FaceAlpha',alpha,'FaceLighting','gouraud');
    % Plot text
    if ~isempty(text_str)
        p = 0.5*(p1+p2);
        h{fig_idx,subfig_idx}.text = text(p(1),p(2),p(3),[' ',text_str],...
            'FontSize',dragzoom_rate*text_fs,'FontName',text_fn,'Color',text_color,'Interpreter',interpreter);
    end
else
    h{fig_idx,subfig_idx}.arrow.Faces = fv.faces;
    h{fig_idx,subfig_idx}.arrow.Vertices = fv.vertices;
    
    if ~isempty(text_str)
        p = 0.5*(p1+p2);
        h{fig_idx,subfig_idx}.text.Position = p;
        h{fig_idx,subfig_idx}.text.String = [' ',text_str];
        h{fig_idx,subfig_idx}.text.FontSize = dragzoom_rate*text_fs;
    end
end
