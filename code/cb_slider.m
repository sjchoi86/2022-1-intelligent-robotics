function cb_slider(~,eventdata,slider_idx)
global g_cb g_edit

% Get the slider value 
slider_val = get(eventdata.AffectedObject, 'Value');
if ~isempty(slider_val)
    g_cb.slider_val{slider_idx} = slider_val;
end

% Write the slider value to the corresponding textbox. 
g_edit{slider_idx}.String = sprintf('%.2f deg',slider_val);

% fprintf('slider_idx:[%d] slider_val:[%.3f].\n',slider_idx,slider_val);
