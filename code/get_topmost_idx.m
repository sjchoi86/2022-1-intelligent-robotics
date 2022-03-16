function idx = get_topmost_idx(chain)
%
% Get the topmost index of a node 
%

idx = round(length(chain.joint)/2); % start from the middle
parent = chain.joint(idx).parent;
while ~isempty(parent)
    idx = parent;
    parent = chain.joint(idx).parent;
    if parent == 0
        break;
    end
end
