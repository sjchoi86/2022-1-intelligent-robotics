function ca
%
% close all with clearning all persistent variables in 'plot_*"
%
close all;
clear_persistents;
end

function clear_persistents
%
% Clear all persistent variables in functions start with 'plot_' (e.g., 'plot_chain')
%
plot_func_names = dir('../**/plot_*.m');
for i_idx = 1:length(plot_func_names)
    plot_func_name = plot_func_names(i_idx).name;
    clear(plot_func_name(1:end-2));
end
end