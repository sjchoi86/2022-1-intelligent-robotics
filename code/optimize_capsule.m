function cap_opt = optimize_capsule(fv)
%
% Optimize capsule 
%

% Reduce once more
reduce_patch_ratio = min(1,600/size(fv.vertices,1)); % #vertex is 100
fv = reducepatch(fv,reduce_patch_ratio);
% Initial guess
radius = max(max(fv.vertices)-min(fv.vertices));
height = max(max(fv.vertices)-min(fv.vertices));
% Optimize
res = 30;
pdata = fv.vertices;
mean_fv = mean(fv.vertices);
use_vector_operation = 1;
x0 = [mean_fv(1) mean_fv(2) -1.5*height ...
    mean_fv(1) mean_fv(2) 1.5*height ...
    radius*2]'; % initial guess
algorithm = 'interior-point';  % 'sqp' 'interior-point'
options = optimoptions('fmincon','Algorithm',algorithm,'display','off',...
    'StepTolerance',1e-10,'OptimalityTolerance',1e-10);
A = []; B = []; Aeq = []; Beq = []; LB = []; UB = [];
NONLCON = @(x)mycon(x,pdata,use_vector_operation);
[x,fval] = fmincon(@(x)myfun(x),x0,A,B,Aeq,Beq,LB,UB,NONLCON,options);
e1 = x(1:3); e2 = x(4:6); r = x(7);
height_opt = norm(e2-e1);
radius_opt = r;
p_opt = (e1 + e2)/2;
diff_vec = e2-e1; diff_unit = diff_vec / norm(diff_vec,'fro');
R_opt = vecRotMat([0 0 1], diff_unit');
cap_opt = get_capsule_shape(...
    'T_offset',pr2t(p_opt,R_opt),'radius',radius_opt,'height',height_opt,'res',res);
cap_opt.p = p_opt;
cap_opt.R = R_opt;
cap_opt.height = height_opt;
cap_opt.radius = radius_opt;
cap_opt.fval = fval;
end

function f = myfun(x)
e1 = x(1:3);
e2 = x(4:6);
r = x(7);
f = norm(e2-e1)*pi*r^2 + 4/3*pi*r^3;
end

function [c,ceq] = mycon(x,pdata,use_vector_operation)
e1 = x(1:3);
e2 = x(4:6);
r = x(7);
Npts = size(pdata,1);
c = zeros(Npts,1);
% use_vector_operation = 1;
if (use_vector_operation)
    % 1.2: tolerance
    tol = 0.95;
    c = dist2segBatch(pdata',e1,e2) - tol*r; % c <= 0
else
    for i=1:Npts
        p = pdata(i,:)';
        c(i) = dist2seg(p,e1,e2) - r; % c <= 0
    end
end
ceq = [];
end

function d = dist2segBatch(pdata,v1,v2)
%DIST2SEG Compute the distance from a point to a line segment
%   Input
%       pdata: an arbitrary point (or NxM point array)
%       v1: start vertex of the segment (Mx1)
%       v2: end vertex of the segment (Mx1)
%   Output
%       d: distance array from the points to the line segment (Nx1)
%
%   http://geomalgorithms.com/a02-_lines.html
v = v2 - v1; % a line segment vector
w = pdata - v1; % a line from starting point to the arbitrary point
c1 = v'*w; % dot(v,w)
c2 = v'*v; % dot(v,v)
Npt = size(pdata,2);
d = zeros(Npt,1);
% the normal vector from the point falls outside the segment
% case 1) closer to v1
ind1 = find(c1<=0);
if(~isempty(ind1))
    d(ind1) = vecnorm(pdata(:,ind1)-v1);
end
% case 2) closer to v2
ind2 = find(c2 <= c1);
if(~isempty(ind2))
    d(ind2) = vecnorm(pdata(:,ind2)-v2);
end
% case 3) the normal vector from the point falls inside the segment
ind3 = setdiff(1:Npt, union(ind1,ind2));
b = c1(ind3)/c2;
pb = v1 + b.*v;

if(~isempty(ind3))
    d(ind3) = vecnorm(pdata(:,ind3)-pb);
end
end


function R = vecRotMat(f,t)
%% Purpose:
%Commonly, it is desired to have a rotation matrix which will rotate one
%unit vector, f,  into another unit vector, t. It is desired to
%find R(f,t) such that R(f,t)*f = t.
%
%This program, vecRotMat is the most
%efficient way to accomplish this task. It uses no square roots or
%trigonometric functions as they are very computationally expensive.
%It is derived from the work performed by Moller and Hughes, which have
%suggested that this method is the faster than any previous transformation
%matrix methods tested.
%
%
%% Inputs:
%f                      [N x 3]                         N number of vectors
%                                                       in which to
%                                                       transform into
%                                                       vector t.
%
%t                      [N x 3]                         N number of vectors
%                                                       in which it is
%                                                       desired to rotate
%                                                       f.
%
%
%% Outputs:
%R                      [3 x 3 x N]                     N number of
%                                                       rotation matrices
%
%% Source:
% Moller,T. Hughes, F. "Efficiently Building a Matrix to Rotate One
% Vector to Another", 1999. http://www.acm.org/jgt/papers/MollerHughes99
%
%% Created By:
% Darin C. Koblick (C) 07/17/2012
% Darin C. Koblick     04/22/2014       Updated when lines are close to
%                                       parallel by checking
%% ---------------------- Begin Code Sequence -----------------------------
%It is assumed that both inputs are in vector format N x 3
dim3 = 2;
%Declare function handles for multi-dim operations
normMD = @(x,y) sqrt(sum(x.^2,y));
anyMD  = @(x) any(x(:));
% Inputs Need to be in Unit Vector Format
if anyMD(single(normMD(f,dim3)) ~= single(1)) || anyMD(single(normMD(t,dim3)) ~= single(1))
    error('Input Vectors Must Be Unit Vectors');
end
%Pre-Allocate the 3-D transformation matrix
R = NaN(3,3,size(f,1));

v = permute(cross(f,t,dim3),[3 2 1]);
c = permute(dot(f,t,dim3),[3 2 1]);
h = (1-c)./dot(v,v,dim3);

idx  = abs(c) > 1-1e-13;
%If f and t are not parallel, use the following computation
if any(~idx)
    %For any vector u, the rotation matrix is found from:
    R(:,:,~idx) = ...
        [c(:,:,~idx) + h(:,:,~idx).*v(:,1,~idx).^2,h(:,:,~idx).*v(:,1,~idx).*v(:,2,~idx)-v(:,3,~idx),h(:,:,~idx).*v(:,1,~idx).*v(:,3,~idx)+v(:,2,~idx); ...
        h(:,:,~idx).*v(:,1,~idx).*v(:,2,~idx)+v(:,3,~idx),c(:,:,~idx)+h(:,:,~idx).*v(:,2,~idx).^2,h(:,:,~idx).*v(:,2,~idx).*v(:,3,~idx)-v(:,1,~idx); ...
        h(:,:,~idx).*v(:,1,~idx).*v(:,3,~idx)-v(:,2,~idx),h(:,:,~idx).*v(:,2,~idx).*v(:,3,~idx)+v(:,1,~idx),c(:,:,~idx)+h(:,:,~idx).*v(:,3,~idx).^2];
end
%If f and t are close to parallel, use the following computation
if any(idx)
    f = permute(f,[3 2 1]);
    t = permute(t,[3 2 1]);
    p = zeros(size(f));
    iidx = abs(f(:,1,:)) <= abs(f(:,2,:)) & abs(f(:,1,:)) < abs(f(:,3,:));
    if any(iidx & idx)
        p(:,1,iidx & idx) = 1;
    end
    iidx = abs(f(:,2,:)) < abs(f(:,1,:)) & abs(f(:,2,:)) <= abs(f(:,3,:));
    if any(iidx & idx)
        p(:,2,iidx & idx) = 1;
    end
    iidx = abs(f(:,3,:)) <= abs(f(:,1,:)) & abs(f(:,3,:)) < abs(f(:,2,:));
    if any(iidx & idx)
        p(:,3,iidx & idx) = 1;
    end
    u = p(:,:,idx)-f(:,:,idx);
    v = p(:,:,idx)-t(:,:,idx);
    rt1 = -2./dot(u,u,dim3);
    rt2 = -2./dot(v,v,dim3);
    rt3 = 4.*dot(u,v,dim3)./(dot(u,u,dim3).*dot(v,v,dim3));
    R11 = 1 + rt1.*u(:,1,:).*u(:,1,:)+rt2.*v(:,1,:).*v(:,1,:)+rt3.*v(:,1,:).*u(:,1,:);
    R12 = rt1.*u(:,1,:).*u(:,2,:)+rt2.*v(:,1,:).*v(:,2,:)+rt3.*v(:,1,:).*u(:,2,:);
    R13 = rt1.*u(:,1,:).*u(:,3,:)+rt2.*v(:,1,:).*v(:,3,:)+rt3.*v(:,1,:).*u(:,3,:);
    R21 = rt1.*u(:,2,:).*u(:,1,:)+rt2.*v(:,2,:).*v(:,1,:)+rt3.*v(:,2,:).*u(:,1,:);
    R22 = 1 + rt1.*u(:,2,:).*u(:,2,:)+rt2.*v(:,2,:).*v(:,2,:)+rt3.*v(:,2,:).*u(:,2,:);
    R23 = rt1.*u(:,2,:).*u(:,3,:)+rt2.*v(:,2,:).*v(:,3,:)+rt3.*v(:,2,:).*u(:,3,:);
    R31 = rt1.*u(:,3,:).*u(:,1,:)+rt2.*v(:,3,:).*v(:,1,:)+rt3.*v(:,3,:).*u(:,1,:);
    R32 = rt1.*u(:,3,:).*u(:,2,:)+rt2.*v(:,3,:).*v(:,2,:)+rt3.*v(:,3,:).*u(:,2,:);
    R33 = 1 + rt1.*u(:,3,:).*u(:,3,:)+rt2.*v(:,3,:).*v(:,3,:)+rt3.*v(:,3,:).*u(:,3,:);
    R(:,:,idx) = [R11 R12 R13; R21 R22 R23; R31 R32 R33];
end
end
