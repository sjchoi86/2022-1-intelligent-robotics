function [m,I_bar] = get_capsule_mass_inertia(cap,varargin)
%
% Get mass and inertia tensor of a capuse
%

% Parse options
ps = inputParser;
addParameter(ps,'density',500);    % density (e.g., 500 [kg/m^3] for firm wood)
parse(ps,varargin{:});
density    = ps.Results.density;

% Compute mass
r = cap.radius;
h = cap.height;
m = density*((4/3)*pi*(r^3)+pi*(r^2)*h); % [kg]

% Get capsule inertia
ixx = m*h*(1/8)*(3*r+2*h);
izz = ixx;
iyy = (2/5)*m*r*r;
I_bar = diag([ixx,iyy,izz]);
