function f = plot_traj(traj,varargin)
%
% Plot a trajectory
%
global zoomPct % zoom rate from dragzoom 
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:10, for j = 1:100, h{i,j}.first_flag = true; end; end; end

% Parse options
ps = inputParser;
addParameter(ps,'fig_idx',1); 
addParameter(ps,'subfig_idx',1);
addParameter(ps,'tlc','k');
addParameter(ps,'tlw',2);
addParameter(ps,'tls','-');
addParameter(ps,'USE_ZOOMRATE',0);
parse(ps,varargin{:});
fig_idx         = ps.Results.fig_idx;
subfig_idx      = ps.Results.subfig_idx;
tlc             = ps.Results.tlc;
tlw             = ps.Results.tlw;
tls             = ps.Results.tls;
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

if isempty(traj)
    traj = nan*[1,1,1];
end

% Plot start here
if h{fig_idx,subfig_idx}.first_flag || ~ishandle(h{fig_idx,subfig_idx}.fig)
    h{fig_idx,subfig_idx}.first_flag = false;
    h{fig_idx,subfig_idx}.fig = figure(fig_idx);
    h{fig_idx,subfig_idx}.traj = plot3(...
        traj(:,1),traj(:,2),traj(:,3),...
        'Color',tlc,'LineWidth',dragzoom_rate*tlw,'LineStyle',tls);
else
    h{fig_idx,subfig_idx}.traj.XData = traj(:,1);
    h{fig_idx,subfig_idx}.traj.YData = traj(:,2);
    h{fig_idx,subfig_idx}.traj.ZData = traj(:,3);
    h{fig_idx,subfig_idx}.traj.LineWidth = dragzoom_rate*tlw;
end

f = h{fig_idx,subfig_idx}.traj;