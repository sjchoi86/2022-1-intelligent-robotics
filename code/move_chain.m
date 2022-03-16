function chain = move_chain(chain,p_root)
%
% Move mocap skeleton
%

root_idx = get_topmost_idx(chain);
chain.joint(root_idx).p = cv(p_root);
chain = fk_chain(chain,'');
