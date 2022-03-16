function plot_title_with_text(title_str,varargin)
%
% Plot a title
%
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:10, h{i}.first_flag = true; end; end

% Parse options
ps = inputParser;
addParameter(ps,'fig_idx',1);
addParameter(ps,'tfs',13);               % title font size
addParameter(ps,'tfn','consolas');       % title font name 
addParameter(ps,'tfc','k');              % title font color
addParameter(ps,'interpreter','none');   % 'none'
parse(ps,varargin{:});
fig_idx         = ps.Results.fig_idx;
tfs             = ps.Results.tfs;
tfn             = ps.Results.tfn;
tfc             = ps.Results.tfc;
interpreter     = ps.Results.interpreter;

% Plot title 
if isequal(interpreter,'latex')
    title_str =  strrep(title_str,'_','-'); % change '_' to '-' for latex formatting
end
if h{fig_idx}.first_flag || ~ishandle(h{fig_idx}.fig)
    h{fig_idx}.first_flag = false;
    h{fig_idx}.fig = figure(fig_idx);
    x = 0.5;
    y = 0.95;
    z = 0.5;
    h{fig_idx}.title = text(x,y,z,title_str,'sc',...
        'HorizontalAlignment','center',...
        'fontsize',tfs,'fontname',tfn,'interpreter',interpreter,'color',tfc);
else
    h{fig_idx}.title.String = title_str;
    h{fig_idx}.title.Color  = tfc;
end
