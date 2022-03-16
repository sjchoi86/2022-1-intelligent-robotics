function idx = idx_cell(cell_data,query_data)
%
% Find the index whose corresponding item matches the 'query_data'
%  We assume that idx is a single or empty index. 
%

idx = '';
for i_idx = 1:length(cell_data)
    if isequal(cell_data{i_idx},query_data)
        idx = i_idx;
        return;
    end
end
