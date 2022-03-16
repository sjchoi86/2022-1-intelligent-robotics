function chain = update_chain_mass_inertia_com(chain,varargin)
%
% Update mass, inertia, and com of link capsules
%

% Parse options
ps = inputParser;
addParameter(ps,'density',500);    % density (e.g., 500 [kg/m^3] for firm wood)
parse(ps,varargin{:});
density    = ps.Results.density;

% Loop
for i_idx = 1:chain.n_link % for all links
    cap_i = chain.link(i_idx).capsule;
    if ~isempty(cap_i)
        [m_i,I_bar_i] = get_capsule_mass_inertia(cap_i,'density',density);
        if ~isfield(chain.link(i_idx),'m')
            chain.link(i_idx).m       = m_i;                    % mass
            chain.link(i_idx).I_bar   = I_bar_i;                % inertia tensor
            chain.link(i_idx).com_bar = t2p(cap_i.T_offset);    % com w.r.t. local coordinates
        elseif chain.link(i_idx).m == 0
            chain.link(i_idx).m       = m_i;                    % mass
            chain.link(i_idx).I_bar   = I_bar_i;                % inertia tensor
            chain.link(i_idx).com_bar = t2p(cap_i.T_offset);    % com w.r.t. local coordinates
        else
            % If the link information is already there, do not update 
        end
    else
        chain.link(i_idx).m       = 0;                      % mass
        chain.link(i_idx).I_bar   = eye(3,3);               % inertia tensor
        chain.link(i_idx).com_bar = cv([0,0,0]);            % com w.r.t. local coordinates
    end
end % for i_idx = 1:chain.n_link % for all links
