function idxs = idxs_cell(cell_data,query_data)
%
% Find the indices whose corresponding items matche the 'query_data'
%

n_query = length(query_data);
idxs = zeros(n_query,1);

for i_idx = 1:n_query
    idxs(i_idx) = idx_cell(cell_data,query_data{i_idx});
end