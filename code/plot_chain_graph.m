function fig = plot_chain_graph(chain,varargin)
%
% Plot kinematic chain graph
%

% Parse input options
ps = inputParser;
addParameter(ps,'fig_idx',2);
addParameter(ps,'fig_pos',[0.5,0.1,0.5,0.3]);
addParameter(ps,'ms',9); % marker size
addParameter(ps,'lw',1); % line width
addParameter(ps,'text_fs',9);
addParameter(ps,'max_str_len',inf);
addParameter(ps,'title_str',chain.name);
addParameter(ps,'title_fs',15);
addParameter(ps,'interpreter','none');
addParameter(ps,'NO_MARGIN',0);
parse(ps,varargin{:});
fig_idx     = ps.Results.fig_idx;
fig_pos     = ps.Results.fig_pos;
ms          = ps.Results.ms;
lw          = ps.Results.lw;
text_fs     = ps.Results.text_fs;
max_str_len = ps.Results.max_str_len;
title_str   = ps.Results.title_str;
title_fs    = ps.Results.title_fs;
interpreter = ps.Results.interpreter;
NO_MARGIN   = ps.Results.NO_MARGIN;

n = chain.n_joint; % number of joint
parents = zeros(1,n);
for i_idx = 1:n
    if ~isempty(chain.joint(i_idx).parent)
        parents(i_idx) = chain.joint(i_idx).parent;
    end
end
[x,y,~] = treelayout(parents); % get tree structure
f = find(parents~=0);
pp = parents(f);
X = [x(f); x(pp); NaN(size(f))];
Y = [y(f); y(pp); NaN(size(f))];
X = X(:);
Y = Y(:);

% Set figure
fig = set_fig(figure(fig_idx),'pos',fig_pos,'AXIS_EQUAL',0,'NO_MARGIN',NO_MARGIN);
x_offset = 0.0;
plot (X+x_offset, Y, 'k-','MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor','w');
for i_idx = 1:n % for all joints
    % Plot node
    if isfield(chain.joint(i_idx),'a') % if rotational axis exists
        if sum(abs(chain.joint(i_idx).a)) == 0
            color = 'w';
        else
            [~,max_idx] = max(abs(chain.joint(i_idx).a'));
            rgb = {'r','g','b'};
            color = rgb{max_idx};
        end
    else
        color = 'w';
    end
    plot (x(i_idx)+x_offset,y(i_idx),'ko',...
        'MarkerSize',ms,'LineWidth',lw,'MarkerFaceColor',color);
    switch interpreter
        case 'latex'
            str = sprintf('~[%02d]%s',i_idx,chain.joint(i_idx).name);
            str = strrep(str,'_','-'); % replace '_' with '-'
        otherwise
            str = sprintf('  [%02d]%s',i_idx,chain.joint(i_idx).name);
    end
    
    % Trim string length
    if length(str) > max_str_len
        str = [str(1:max_str_len),'..'];
    end
    
    % Text joint name
    text(x(i_idx)+x_offset,y(i_idx),str,'FontSize',text_fs,...
        'interpreter',interpreter,'FontName','Consolas');
end

% Title
switch interpreter
    case 'latex'
        title_str = strrep(title_str,'_','-'); % replace '_' with '-'
end
title(title_str,'fontsize',title_fs,'fontname','Consolas','interpreter',interpreter);
% plot_title_with_text(title_str,'fig_idx',fig_idx,...
%     'tfs',title_fs,'tfn','Consolas','interpreter',interpreter);
xlim([0,1]);
ylim([0,1]);
% xlim([min(x)-x_margin,1+x_margin]);
% ylim([-y_margin/2,1+y_margin]);
set(gcf,'Color','w'); % background color to white 
axis off;
dragzoom;
