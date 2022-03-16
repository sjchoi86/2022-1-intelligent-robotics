function set_sliders_with_text_boxes(fig,slider_names,slider_values,varargin)
%
% Set sliders with text boxes
%
% Slider values can be read by
% for i_idx = 1:length(g_cb.slider_val)
%     val_i = g_cb.slider_val{i_idx};
% end
global g_cb g_edit g_slider

% Parse options
ps = inputParser;
addParameter(ps,'sliderstep',1/30); 
addParameter(ps,'slider_min',-180); 
addParameter(ps,'slider_max',180);
addParameter(ps,'bg_colors','');
parse(ps,varargin{:});
sliderstep = ps.Results.sliderstep;
slider_min = ps.Results.slider_min;
slider_max = ps.Results.slider_max;
bg_colors  = ps.Results.bg_colors;

n_slider = length(slider_names); % number of sliders
for i_idx = 1:n_slider
    joint_name = slider_names{i_idx};
    joint_position = slider_values(i_idx,1);
    
    % Common param
    pos_ymin = 0.03+0.92*(1/n_slider)*(n_slider-i_idx);
    edit_height = 1/n_slider;
    
    % Add a text box for showing joint name
    pos_xmin = 0.05;
    edit_len = 0.14;
    if isempty(bg_colors)
        bg_color = [0.5,0.7,0.9,0.5];
    else
        bg_color = bg_colors(i_idx,:);
    end
    h_edit{i_idx} = uicontrol('parent',fig,...
        'style','edit','units','normalized',...
        'position',[pos_xmin,pos_ymin,edit_len,edit_height],...
        'BackgroundColor', bg_color,...
        'String',sprintf('[%d]%s',i_idx,joint_name) ...
        );
    
    % Add a slidler for changing joint value
    pos_xmin = 0.2;
    slider_len = 0.65;
    slider_height = edit_height;
    g_slider{i_idx} = uicontrol('parent',fig,...
        'style','slider','units','normalized',...
        'position',[pos_xmin,pos_ymin,slider_len,slider_height],...
        'sliderstep',[sliderstep,sliderstep],...
        'BackgroundColor',bg_color,...
        'Value',joint_position,...
        'min',slider_min,'max',slider_max);
    callback = @(a,b) cb_slider(a,b,i_idx);
    addlistener(g_slider{i_idx},'Value','PostSet',callback); % callback listener
    g_cb.slider_val{i_idx} = joint_position;
    
    % Add a text box for showing joint value in slider
    pos_xmin = 0.85;
    edit_len = 0.1;
    g_edit{i_idx} = uicontrol('parent',fig,...
        'style','edit','units','normalized',...
        'position',[pos_xmin,pos_ymin,edit_len,edit_height],...
        'BackgroundColor',bg_color,...
        'String',sprintf('%.2f deg',joint_position)...
        );
    
end
