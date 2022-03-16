function ik_plot_info = get_ik_plot_info_from_ik_info(ik_info)
%
% Get 'ik_plot_info' from 'ik_info'
%
ik_joint_names = cell(1,ik_info.n_trgt);
ik_types = cell(1,ik_info.n_trgt);
ik_trgt_coords = cell(1,ik_info.n_trgt);
for ik_idx = 1:ik_info.n_trgt
    ik_joint_names{ik_idx} = ik_info.trgt_joint_names{ik_idx};
    ik_types{ik_idx} = ik_info.trgt_types{ik_idx};
    ik_trgt_coords{ik_idx} = ik_info.trgt_coords{ik_idx};
end
ik_plot_info.joint_names = ik_joint_names;
ik_plot_info.types = ik_types;
ik_plot_info.trgt_coords = ik_trgt_coords;
