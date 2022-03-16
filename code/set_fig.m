function fig = set_fig(fig,varargin)
%
% Set figure
%
% Usage:
% set_fig(figure(1),'pos',[0.5,0.4,0.3,0.55],...
%     'view_info',view_info,'axis_info',1.5*[-1,+1,-1,+1,-1,+1],'AXIS_EQUAL',1,'GRID_ON',1,...
%     'REMOVE_MENUBAR',1,'USE_DRAGZOOM',1,'SET_CAMLIGHT',1,'SET_MATERIAL','METAL',...
%     'SET_AXISLABEL',1,'afs',18);
%
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:100,h{i}.first_flag = true; end; end

% Parse options
ps = inputParser;
addParameter(ps,'pos',[0.0,0.5,0.3,0.5]);            % position of the figure
addParameter(ps,'view_info','');                     % view information
addParameter(ps,'axis_info','');                     % axis information
addParameter(ps,'HOLD_ON',1);                        % hold on
addParameter(ps,'AXIS_EQUAL',1) ;                    % axis equal
addParameter(ps,'AXIS_OFF',0);                       % axis off
addParameter(ps,'GRID_ON',1);                        % grid on
addParameter(ps,'REMOVE_MENUBAR',1);                 % remove menubar
addParameter(ps,'USE_DRAGZOOM',1);                   % use dragzoome
addParameter(ps,'SET_CAMLIGHT',1);                   % set camera light
addParameter(ps,'SET_MATERIAL','METAL');             % set material property 'SHINY' 'DULL' 'METAL'
addParameter(ps,'SET_AXISLABEL',1);                  % axis label
addParameter(ps,'ax_str','X');                       % x label string
addParameter(ps,'ay_str','Y');                       % y label string
addParameter(ps,'az_str','Z');                       % z label string
addParameter(ps,'afs',13);                           % axis font size
addParameter(ps,'interpreter','latex');               % text interpreter
addParameter(ps,'gcf_color','w');                    % background color
addParameter(ps,'NO_MARGIN',0);                      % no margin on axes
parse(ps,varargin{:});
pos             = ps.Results.pos;
view_info       = ps.Results.view_info;
axis_info       = ps.Results.axis_info;
HOLD_ON         = ps.Results.HOLD_ON;
AXIS_EQUAL      = ps.Results.AXIS_EQUAL;
AXIS_OFF        = ps.Results.AXIS_OFF;
GRID_ON         = ps.Results.GRID_ON;
REMOVE_MENUBAR  = ps.Results.REMOVE_MENUBAR;
USE_DRAGZOOM    = ps.Results.USE_DRAGZOOM;
SET_CAMLIGHT    = ps.Results.SET_CAMLIGHT;
SET_MATERIAL    = ps.Results.SET_MATERIAL;
SET_AXISLABEL   = ps.Results.SET_AXISLABEL;
ax_str          = ps.Results.ax_str;
ay_str          = ps.Results.ay_str;
az_str          = ps.Results.az_str;
afs             = ps.Results.afs;
interpreter     = ps.Results.interpreter;
gcf_color       = ps.Results.gcf_color;
NO_MARGIN       = ps.Results.NO_MARGIN;

% Set figure configuration
if h{fig.Number}.first_flag || ~ishandle(h{fig.Number}.fig)
    h{fig.Number}.first_flag = false;
    h{fig.Number}.fig = fig;
    
    % No margin
    if NO_MARGIN
        axes('Parent',fig,'Position',[0,0,1,1]);
    end
    
    % Set figure position
    if ~isempty(pos)
        sz = get(0, 'ScreenSize');
        fig_pos = [pos(1)*sz(3),pos(2)*sz(4),pos(3)*sz(3),pos(4)*sz(4)];
        set(fig,'Position',fig_pos);
    end
    if HOLD_ON
        hold on;
    end
    
    % View information
    if ~isempty(view_info)
        switch length(view_info)
            case 1, view(view_info);
            case 2, view(view_info(1),view_info(2));
            case 3, view(view_info);
        end
    end
    
    % Make axis equal
    if AXIS_EQUAL
        axis equal;
    end
    
    % Axis information (this should come AFTER axis equal)
    if ~isempty(axis_info)
        axis(axis_info);
    end
    
    % Axis off
    if AXIS_OFF
        axis off;
    end
    
    % Grid on
    if GRID_ON
        grid on;
    end
    
    % Remove menubar
    if REMOVE_MENUBAR
        set(h{fig.Number}.fig,'MenuBar','none');
    end
    
    % Use dragzoom
    if USE_DRAGZOOM
        dragzoom;
    end
    
    % Set camera light
    if SET_CAMLIGHT
        camlight('infinite');
    end
    
    % Set material
    if ~isempty(SET_MATERIAL)
        material(SET_MATERIAL);
    end
    
    % Set axis label
    if SET_AXISLABEL
        xlabel(ax_str,'fontname','consolas','fontsize',afs,'interpreter',interpreter);
        ylabel(ay_str,'fontname','consolas','fontsize',afs,'interpreter',interpreter);
        zlabel(az_str,'fontname','consolas','fontsize',afs,'interpreter',interpreter);
    end
    
    % Set background color
    if ~isempty(gcf_color)
        set(gcf,'color',gcf_color);
    end
    
    
else
    
end

fig = h{fig.Number}.fig;
