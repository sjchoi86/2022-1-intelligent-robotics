function sc_checks = get_sc_checks(chain,varargin)
%
% Get link indices to check self-collision
%

% Parse options
ps = inputParser;
addParameter(ps,'collision_margin',max(chain.sz.xyz_len)/100);
parse(ps,varargin{:});
collision_margin = ps.Results.collision_margin;

sc_checks = [];
for i_idx = 1:chain.n_link % for all links
    link_i = chain.link(i_idx);
    cap_i = link_i.capsule;
    if isempty(cap_i), continue; end % ignore empty capsule
    joint_idx_i = link_i.joint_idx;
    if isempty(joint_idx_i) 
        joint_idx_i = get_topmost_idx(chain);
    end
    T_i = pr2t(chain.joint(joint_idx_i).p,chain.joint(joint_idx_i).R);
    cl_i = get_capsule_line(T_i,cap_i);
    for j_idx = i_idx+1:chain.n_link
        link_j = chain.link(j_idx);
        cap_j = link_j.capsule;
        if isempty(cap_j), continue; end % ignore empty capsule
        joint_idx_j = link_j.joint_idx;
        if isempty(joint_idx_j)
            joint_idx_j = get_topmost_idx(chain);
        end
        T_j = pr2t(chain.joint(joint_idx_j).p,chain.joint(joint_idx_j).R);
        cl_j = get_capsule_line(T_j,cap_j);
        line_dist = get_dist_lines(cl_i.p1,cl_i.p2,cl_j.p1,cl_j.p2);
        cap_dist = line_dist - cap_i.radius - cap_j.radius;
        if cap_dist > collision_margin
            sc_checks = cat(1,sc_checks,[i_idx,j_idx]);
        end
    end
end

