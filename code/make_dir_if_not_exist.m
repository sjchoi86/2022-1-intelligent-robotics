function make_dir_if_not_exist(folder_path)
%
% Make a directory if not exist
%
if ~exist(folder_path,'dir')
    mkdir(folder_path);
    fprintf(2,'mkdir [%s].\n',folder_path);
end
