function [SC,sc_link_pairs] = check_sc(chain,varargin)
%
% Check self-collision
%

% Parse options
ps = inputParser;
addParameter(ps,'collision_margin',max(chain.sz.xyz_len)/100);
parse(ps,varargin{:});
collision_margin = ps.Results.collision_margin;

% Check self-collision
SC = 0; sc_link_pairs = [];
for sc_idx = 1:size(chain.sc_checks,1)
    sc_check = chain.sc_checks(sc_idx,:);
    i_idx = sc_check(1); j_idx = sc_check(2);
    % Link i
    link_i = chain.link(i_idx);
    cap_i = link_i.capsule;
    joint_idx_i = link_i.joint_idx;
    T_i = pr2t(chain.joint(joint_idx_i).p,chain.joint(joint_idx_i).R);
    cl_i = get_capsule_line(T_i,cap_i);
    % Link j
    link_j = chain.link(j_idx);
    cap_j = link_j.capsule;
    joint_idx_j = link_j.joint_idx;
    T_j = pr2t(chain.joint(joint_idx_j).p,chain.joint(joint_idx_j).R);
    cl_j = get_capsule_line(T_j,cap_j);
    % Compute the distance between two capsules
    line_dist = get_dist_lines(cl_i.p1,cl_i.p2,cl_j.p1,cl_j.p2);
    cap_dist = line_dist - cap_i.radius - cap_j.radius;
    if cap_dist < -collision_margin
        SC = 1;
        sc_link_pairs = cat(1,sc_link_pairs,sc_check);
    end
end
