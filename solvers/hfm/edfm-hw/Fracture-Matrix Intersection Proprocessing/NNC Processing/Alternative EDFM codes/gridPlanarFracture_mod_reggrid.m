function Gf = gridPlanarFracture_mod1(G_matrix, fracplane, tol, varargin)
% Modifications: calculation for theta was changed to theta = -atan(diffp(1,2)/diffp(1,1));
% Removed lines 252-254. Changed the input to be explicitly G_matrix.
% Verified that fracture aperture is used correctly. Have not removed the
% Cartesian & PEBI portions (lazy). Fixed triangle meshing part. If not
% rectangular, the last argument needs to have the first vertex repeated,
% i.e. [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, xyp, [xyp;xyp(1,:)]);
% Remove cartesian and pebi sections.
%
% This function can be used to grid a planar fracture using a cartesian,
% triangular or PEBI mesh. See preProcessingFractures for input options.
%
% SYNOPSIS:
%   Gf = gridPlanarFracture(G, fracplane, scaledplane)
%   Gf = gridPlanarFracture(G, fracplane, scaledplane, 'pn1', pv1)
%
% REQUIRED PARAMETERS:
%
%   G_matrix    - Matrix grid data structure containing geometrical information.
%
%   fracplane   - Structure containing a set of coplanar points that define
%                 the fracture polygon and a value for the fracture
%                 aperture.
%
%   scaledplane - Same as fracplane with geometrical information in
%                 dimensionless form.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   
%   cellSize           - Dimensionless element size (>0 and <1) for the
%                        fracture grid.
%                
%   minTriangles       - Can be used to set a minimum on the number of
%                        cells for a triangular grid.
%
% RETURNS:
%   Gf - Fracture grid structure with geometrical information.

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('cellSize'     ,   -1        , ...
             'minTriangles' ,   10       );
opt = merge_options(opt, varargin{:});

% Rotate points about the origin to align with the xy plane
xyp = rotatePlane(fracplane.points,[0 0 1]);

% Extract xy-coordinates and store the z coordinate
z = xyp(1,3);
xyp = xyp(:,1:2);

% % Define a scale. Save the scale. Compute element size
% maxpts=max(xyp);
% minpts=min(xyp);
% scales=maxpts-minpts;
% scale=min(scales);
% xyp=xyp/scale;

% Rotate axis to align 1 edge of the planar polygon with the x-axis 
diffp = diff(xyp);
theta = -atan(diffp(1,2)/diffp(1,1));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
xyp = transpose(R*xyp');

% Default scaled cellsize=0.25
if opt.cellSize < tol
    h = 0.25;
else
    h = opt.cellSize;
end

% Refine h or average element size
% h = (h<1)*h^2 + (h>1)*sqrt(h); % assuming esize<0
% hstr = num2str(h,'%f');
% sigfigs = hstr(hstr~='0' & hstr~='.'); % ex: if h = 1.43, sigfigs = 143 
% sigfigs2 = find(hstr~='0' & hstr~='.')-1; % ex: if h = 0.04, sigfigs2 = 3
% nsigfigs = max(length(sigfigs),sigfigs2(1)); 
% if nsigfigs>5, nsigfigs = nsigfigs-2; end
% h = round(h,nsigfigs);
% if h<0.01, h = 0.01; elseif h>0.5, h=0.5; end

% Expand bounding box for triangulation
expand = h*10;
bbox = [min(xyp(:,1)) - expand, min(xyp(:,2)) - expand; ...
    max(xyp(:,1)) + expand, max(xyp(:,2)) + expand];


%-----------------------------Triangle------------------------------------%
count = 0;
iter = 100;
while true
    count = count + 1;
    assert(count<10,' ');
    try
        [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, xyp, [xyp;xyp(1,:)]);
        
        assert(size(t,1)>opt.minTriangles,'');
        close gcf
        Gf = triangleGrid(p,t);
        
        % Rotate axis back to original position
        nc = Gf.nodes.coords;
        R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        nc = transpose(R*nc');
        Gf.nodes.coords = nc;
        
%         computeGeometry(Gf);
%         computeGeometry(makeLayeredGrid(Gf,1));
        break;
    catch
        % Perhaps element size is too big, make it smaller.
        h = h/2;
        iter = iter+50;  
    end
end

%----------------------------Gridding complete----------------------------%

% % Scale back
% Gf.nodes.coords=Gf.nodes.coords*scale;

% Add z-coordinate of the original set of points post rotation
nc = [Gf.nodes.coords,repmat(z,Gf.nodes.num,1)];
if isfield(Gf,'cartDims')
    cdims = Gf.cartDims;
    Gf = computeGeometry(makeLayeredGrid(Gf,1));
    Gf.cartDims = cdims;
else
    Gf = computeGeometry(makeLayeredGrid(Gf,1));
end

unitNormal = fracplane.normal/norm(fracplane.normal);
tdist = fracplane.aperture/2;

% Rotate points about the origin back to their original position
use_points = rotatePlane(nc,+unitNormal,'useNormal',[0 0 1]);

% Expand along normal to create an aperture
pminus = use_points - repmat(unitNormal,size(use_points,1),1)*tdist; %correct
pplus = use_points + repmat(unitNormal,size(use_points,1),1)*tdist; %correct
nc = [pminus;pplus];

% Add new node coordinates to fracture grid and re-compute geometry
Gf.nodes.coords = nc;
Gf = computeGeometry(Gf);
Gf.cells.volumes = abs(Gf.cells.volumes);
Gf.faces.areas = abs(Gf.faces.areas);
return

