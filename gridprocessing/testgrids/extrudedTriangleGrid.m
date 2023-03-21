function [G, varargout] = extrudedTriangleGrid(maxarea, varargin)
%Build a synthetic grid with a curved fault in the lateral direction
%
% SYNOPSIS:
%   G = extrudedTriangleGrid()
%   G = extrudedTriangleGrid(maxarea)
%   G = extrudedTriangleGrid(maxarea, dual)
%
% PARAMETERS:
%   maxarea - maximal area of triangles in the areal grid from which the
%             prismatic grid is extruded
%             Default: 100
%   dual    - if true, use a lateral pebi grid rather than a triangle grid
%             as basis for the extrusion process.
%             Default: false
%
% RETURNS:
%   G       - a standard grid structure
%
% EXAMPLES:
%
%   % Prismatic grid
%   G = extrudedTriangleGrid(true);
%   plotGrid(G); view(-30,60), axis tight off
%
%   % PEBI grid
%   G = extrudedTriangleGrid(50, true);
%   plotGrid(G); view(-30,60), axis tight off

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

mlist = mrstModule();
mrstModule add triangle;

if nargin < 1, maxarea = 100; end
   nlayers = 15;

   pslg  = createOutline();

   if true,
      fault       = createFault();
      pslg.edges  = [pslg.edges;  fault.edges+size(pslg.points,1)];
      pslg.points = [pslg.points; fault.points];
   end

   G = mex_triangleGrid(pslg.points, pslg.edges, 'maxArea', maxarea, ...
                        'segmentmarkers', int32(1:size(pslg.edges, 1)));

   if nargin == 2,
      % Make an areal pebi grid rather than the triangular grid.
      G = pebi(G);
   end

   H = computeGeometry(G);
   n = G.cells.num;

   % Find cells to left/right of fault.
   c = findCells(H, fault);

   if (nargin > 1) && (nargout > 1),
      varargout{1} = findEnclosingCell(H, varargin{1});
   end

   G = makeLayeredGrid(G, nlayers);
   G = transformGrid(G);

   c  = repmat(c, [nlayers, 1]);
   c1 = [repmat(true ([n, 1]), [ 3          , 1]); ...
         repmat(false([n, 1]), [nlayers - 3 , 1])];

   c2 = [repmat(false([n, 1]), [nlayers - 5 , 1]); ...
         repmat(true ([n, 1]), [ 5          , 1])];

   if isfield(G.faces, 'tag')
      G.faces = rmfield(G.faces, 'tag');
   end

   G = removeCells(G, (c & c1) | ((~ c) & c2));
   mrstModule('reset', mlist{:});
end

%--------------------------------------------------------------------------

function G = transformGrid(G)
   x = G.nodes.coords(:,1);
   y = G.nodes.coords(:,2);
   z = G.nodes.coords(:,3);

   fun = @(x,y,t) -0.05*sin(3*t*pi*(x-0.5)) -0.075*sin(pi*(y+2*x));
   n = G.nodes.num/(G.numLayers+1);

   I = (1:n);
   xmax = max(x);ymax=max(y);
   z(I) = 40*fun(x(I)/xmax, y(I)/ymax,rand);

   for i=2:G.numLayers+1,
      I = (i-1)*n+(1:n);
      z(I) = z(I-n) + (0.5+rand);%*fun(x(I), y(I),rand);
   end

   G.nodes.coords = [x,y,z];
end

%--------------------------------------------------------------------------

function c = findCells(G, fault)
   c = G.cells.centroids(:,1) < ppval(fault.ppx, G.cells.centroids(:,2));
end

%--------------------------------------------------------------------------

function outline = createOutline()
   outline.points = [0,0; 1,0; 1,1; 0,1]*100;
   outline.edges  = [1,2; 2,3; 3,4; 4,1];
end

%--------------------------------------------------------------------------

function fault = createFault()
   pts   = [0.5 , 0.0 ; ...
            0.5 , 0.3 ; ...
            0.4 , 0.6 ; ...
            0.4 , 1.0 ].*100;

   fault = splineRefine2(pts);
   n     = size(fault.points, 1);

   fault.edges = [(1:n-1) .', (2:n) .'];
end

%--------------------------------------------------------------------------

function curve = splineRefine2(pts)
   ppx = spline(pts(:,2), pts(:,1));

   Y = linspace(0, 100, 20)';
   X = ppval(ppx, Y);

   curve.points = [X,Y];
   curve.ppx = ppx;
end

%--------------------------------------------------------------------------

%function curve = splineRefine(pts)
%   t   = linspace(0, 1, size(pts,1));
%   ppx = spline(t, pts(:,1));
%   ppy = spline(t, pts(:,2));
%
%   tt   = linspace(0, 1, 20)';
%   X = ppval(ppx, tt);
%   Y = ppval(ppy, tt);
%   curve.points = [X,Y];
%   curve.ppx = ppx;
%   curve.ppy = ppy;
%end
