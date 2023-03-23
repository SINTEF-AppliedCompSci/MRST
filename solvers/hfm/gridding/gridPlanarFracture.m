function Gf = gridPlanarFracture(G, fracplane, scaledplane, varargin)
% This function can be used to grid a planar fracture using a cartesian,
% triangular or PEBI mesh. See preProcessingFractures for input options.
%
% SYNOPSIS:
%   Gf = gridPlanarFracture(G, fracplane, scaledplane)
%   Gf = gridPlanarFracture(G, fracplane, scaledplane, 'pn1', pv1)
%
% REQUIRED PARAMETERS:
%
%   G           - Grid data structure containing geometrical information.
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
%   gridType           - Type of fracture mesh desired. Possible options:
%                        (a) 1 - Cartesian mesh: Possible only when the
%                        fracture plane is rectangular.
%
%                        (b) 2 - Triangle mesh: Uses the external library
%                        'distmesh' to create the nodes and connectivity
%                        list required to create a triangular grid.
%
%                        (c) 3 - PEBI mesh: Voronoi grid. External
%                        library 'distmesh' required to create underlying
%                        triangulation.
%
%   cellSize           - Dimensionless element size (>0 and <1) for the
%                        fracture grid.
%                
%   rectangular        - Indicates if the fracture plane is a rectangle. 
%
%   minTriangles       - Can be used to set a minimum on the number of
%                        cells for a triangular grid.
%
%   scale              - Scaling factor for fracture cellSize. Useful when
%                        physical dimensions vary by orders of magnitude in
%                        each direction. Not used if cellSize is provided
%                        as input.
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

opt = struct('gridType'     ,   2        , ...
             'cellSize'     ,   0        , ...
             'rectangular'  ,   false    , ...
             'minTriangles' ,   10       , ...
             'scale'        ,   2        );
opt = merge_options(opt, varargin{:});

% Rotate points about the origin to align with the xy plane
xyp = rotatePlane(scaledplane.points,[0 0 1]);

% Extract xy-coordinates and store the z coordinate
z = xyp(1,3);
xyp = xyp(:,1:2);

% Rotate axis to align 1 edge of the planar polygon with the x-axis 
diffp = diff(xyp);
theta = abs(atan(diffp(1,2)/diffp(1,1)));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
xyp = transpose(R*xyp');

if opt.cellSize == 0
    h = 1;
    h = h/opt.scale;
else
    h = opt.cellSize;
end

% Refine h or average element size
h = (h<1)*h^2 + (h>1)*sqrt(h); % assuming esize<0
hstr = num2str(h,'%f');
sigfigs = hstr(hstr~='0' & hstr~='.'); % ex: if h = 1.43, sigfigs = 143 
sigfigs2 = find(hstr~='0' & hstr~='.')-1; % ex: if h = 0.04, sigfigs2 = 3
nsigfigs = max(length(sigfigs),sigfigs2(1)); 
if nsigfigs>5, nsigfigs = nsigfigs-2; end
h = round(h,nsigfigs);
if h<0.01, h = 0.01; elseif h>0.5, h=0.5; end


% Expand bounding box for triangulation
expand = h*10;
bbox = [min(xyp(:,1)) - expand, min(xyp(:,2)) - expand; ...
    max(xyp(:,1)) + expand, max(xyp(:,2)) + expand];

%-----------------------------Cartesian-----------------------------------%
count = 0;
iter = 100;
if opt.rectangular && opt.gridType == 1
    lx = max(xyp(:,1)) - min(xyp(:,1));
    ly = max(xyp(:,2)) - min(xyp(:,2));
    nx = ceil(1/h/lx);
    ny = ceil(1/h/ly);
    x = linspace(min(xyp(:,1)),max(xyp(:,1)),nx);
    y = linspace(min(xyp(:,2)),max(xyp(:,2)),ny);
    Gf = tensorGrid(x,y);
    
    % Rotate axis back to original position
    nc = Gf.nodes.coords;
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    nc = transpose(R*nc');
    Gf.nodes.coords = nc;
    
    Gf = computeGeometry(Gf);
elseif ~opt.rectangular && opt.gridType == 1
    fprintf('\nFracture is not rectangular, trying to generate a triangle grid.\n'); 
    opt.gridType = 2;
end
%-----------------------------Triangle------------------------------------%
while true && opt.gridType~=1
    count = count + 1;
    assert(count<10,' ');
    try
        if opt.rectangular
            fd = @ (p) drectangle(p, min(xyp(:,1)), max(xyp(:,1)), min(xyp(:,2)), max(xyp(:,2)));
            [p,t] = distmesh_2d(fd, @huniform, h, bbox, iter, xyp);
        else
            [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, xyp, xyp);
        end
        assert(size(t,1)>opt.minTriangles,'');
        close gcf
        Gf = triangleGrid(p,t);
        
        % Rotate axis back to original position
        nc = Gf.nodes.coords;
        R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        nc = transpose(R*nc');
        Gf.nodes.coords = nc;
        
        computeGeometry(Gf);
        computeGeometry(makeLayeredGrid(Gf,1));
        break;
    catch
        % Perhaps element size is too big, make it smaller.
        h = h/2;
        iter = iter+50;  
    end
end
%---------------------------------PEBI------------------------------------%
iter = 50;
if opt.gridType == 3
    while true
        count = count + 1;
        try
            Gf = pebi(Gf);
            close gcf
            assert(size(t,1)>opt.minTriangles*5,'');         
            computeGeometry(Gf);
            break;
        catch
            if count>=4
                opt.gridType = 2;
                fprintf('\nPEBI grid could not be generated successfully. Trying to generate a triangle grid.\n');
                break;
            end
            h = h/2;
            iter = iter+50;
            if opt.rectangular
                fd = @ (p) drectangle(p, min(xyp(:,1)), max(xyp(:,1)), min(xyp(:,2)), max(xyp(:,2)));
                [p,t] = distmesh_2d(fd, @huniform, h, bbox, iter, xyp);
            else
                [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, xyp, xyp);
            end
            Gf = triangleGrid(p,t);
        end
    end
end
%--------------------Triangle if PEBI fails repeatedly--------------------%
if opt.gridType == 2 && count>=4
    iter = 50; h = h*2^(count); count = 0; 
    while true
        count = count + 1;
        assert(count<10,' ');
        try
            if opt.rectangular
                fd = @ (p) drectangle(p, min(xyp(:,1)), max(xyp(:,1)), min(xyp(:,2)), max(xyp(:,2)));
                [p,t] = distmesh_2d(fd, @huniform, h, bbox, iter, xyp);
            else
                [p,t] = distmesh_2d(@dpoly, @huniform, h, bbox, iter, xyp, xyp);
            end
            close gcf
            assert(size(t,1)>opt.minTriangles,'');
            Gf = triangleGrid(p,t);
            
            % Rotate axis back to original position
            nc = Gf.nodes.coords;
            R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            nc = transpose(R*nc');
            Gf.nodes.coords = nc;
            
            computeGeometry(Gf);
            computeGeometry(makeLayeredGrid(Gf,1));
            break;
        catch
            % Perhaps element size is too big, make it smaller.
            h = h/2;
            iter = iter+50;
            
        end
    end
end
%----------------------------Gridding complete----------------------------%

% Add z-coordinate of the original set of points post rotation
nc = [Gf.nodes.coords,repmat(z,Gf.nodes.num,1)];
if isfield(Gf,'cartDims')
    cdims = Gf.cartDims;
    Gf = computeGeometry(makeLayeredGrid(Gf,1));
    Gf.cartDims = cdims;
else
    Gf = computeGeometry(makeLayeredGrid(Gf,1));
end


unitNormal = scaledplane.normal/norm(scaledplane.normal);
tdist = fracplane.aperture/2;

% Rotate points about the origin back to their original position
use_points = rotatePlane(nc,unitNormal,'useNormal',[0 0 1]);

% Scale back to dimensional space
dims = max(G.nodes.coords);
use_points = [use_points(:,1)*dims(1) use_points(:,2)*dims(2) use_points(:,3)*dims(3)];

% Expand along normal to create an aperture
pminus = use_points - repmat(unitNormal,size(use_points,1),1)*tdist;
pplus = use_points + repmat(unitNormal,size(use_points,1),1)*tdist;
nc = [pminus;pplus];

% Add new node coordinates to fracture grid and re-compute geometry
Gf.nodes.coords = nc;
Gf = computeGeometry(Gf);
Gf.cells.volumes = abs(Gf.cells.volumes);
Gf.faces.areas = abs(Gf.faces.areas);
return
