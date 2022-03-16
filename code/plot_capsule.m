function plot_capsule(cap,varargin)
%
% Plot a capsule
%
persistent h

% Make enough handlers at the first
if isempty(h), for i = 1:10, for j = 1:100, h{i,j}.first_flag = true; end; end; end

% Parse input arguments
ps = inputParser;
addParameter(ps,'fig_idx',1);
addParameter(ps,'subfig_idx',1);
addParameter(ps,'T',pr2t(cv([0,0,0]),eye(3,3)));
addParameter(ps,'cfc','r');              % capsule face color
addParameter(ps,'cfa',0.5);              % capsule face alpha
addParameter(ps,'cec','k');              % capsule edge color
addParameter(ps,'cea',0.5);              % capsule edge alpha
addParameter(ps,'RESET',0);
parse(ps,varargin{:});
T               = ps.Results.T;
fig_idx         = ps.Results.fig_idx;
subfig_idx      = ps.Results.subfig_idx;
cfc             = ps.Results.cfc;
cfa             = ps.Results.cfa;
cec             = ps.Results.cec;
cea             = ps.Results.cea;
RESET           = ps.Results.RESET;

if RESET % reset capsule
    for fig_idx = 1:10
        for subfig_idx = 1:100
            if ~h{fig_idx,subfig_idx}.first_flag
                delete(h{fig_idx,subfig_idx}.cap_patch); % delete capsue
                h{fig_idx,subfig_idx}.first_flag = true; % revert flag
            end
        end
    end
    return;
end

if h{fig_idx,subfig_idx}.first_flag || (~ishandle(h{fig_idx,subfig_idx}.fig))
    h{fig_idx,subfig_idx}.first_flag = false;
    h{fig_idx,subfig_idx}.fig = figure(fig_idx);
    % Plot a capsule
    h{fig_idx,subfig_idx}.cap_patch = ...
        patch('faces',cap.faces,'vertices',cap.vertices,...
        'FaceColor',cfc,'EdgeColor',cec,'EdgeAlpha',cea,...
        'FaceLighting','gouraud','AmbientStrength', 0.5, 'FaceAlpha', cfa);
    h{fig_idx,subfig_idx}.cap_patch_t = hgtransform;
    set(h{fig_idx,subfig_idx}.cap_patch,...
        'parent',h{fig_idx,subfig_idx}.cap_patch_t);
    tform = T*cap.T_offset;
    set(h{fig_idx,subfig_idx}.cap_patch_t,'Matrix',tform);
else
    tform = T*cap.T_offset;
    if isnan(sum(tform,'all'))
        delete(h{fig_idx,subfig_idx}.cap_patch);
    else
        set(h{fig_idx,subfig_idx}.cap_patch_t,'Matrix',tform);
    end
end


