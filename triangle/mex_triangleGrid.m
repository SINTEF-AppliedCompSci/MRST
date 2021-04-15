function G = mex_triangleGrid(points, edges, varargin)
%Construct 2d triangle grid in physical space.
%
% SYNOPSIS:
%   G = mex_triangleGrid(pointlist, edgelist)
%   G = mex_triangleGrid(pointlist, edgelist, 'pn1', pv1, ...)
%
% PARAMETERS:
%   pointlist -
%            List of vertex coordinates.  Must be an m-by-2 (double) array
%            of (X,Y) coordinate tuples--one tuple for each vertex.
%
%   edgelist -
%            List of edges defining the boundary of the domain.  Must be an
%            n-by-2 integer array of start- and end nodes (vertices) of
%            individual edges.
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.
%            The supported options are:
%              maxArea  -- Maximum area (m^2) of individual triangles.
%
%              minAngle -- Minimum angle in triangles (degrees).
%                          Use sensible values here (0 < minAngle < 40)
%                          lest the Triangle software fail to compute a
%                          triangulation.
%
%              verbose  -- Whether or not to display progress information
%                          while triangulating the domain.
%                          Logical.  Default value: verbose = false.
%
% RETURNS:
%   G - Grid structure as detailed in grid_structure, without geometric
%       primitives.  Use function computeGeometry to compute those values.
%
% NOTE:
%   This function invokes the Triangle software package.   See website
%
%       http://www.cs.cmu.edu/~quake/triangle.html
%
%   for availability and terms and conditions for use.
%
% EXAMPLE:
%   % Make a 10m-by-5m grid.
%   points = [ 0 , 0 ; 5 , 0 ; 5 , 10 ; 0 , 10 ];
%   edges  = [ 1 , 2 ; 2 , 3 ; 3 ,  4 ; 4 ,  1 ];
%   G = mex_triangleGrid(points, edges, 'maxArea', 0.3);
%
%   % Plot the grid in 3D-view.
%   f = plotGrid(G); axis equal tight; view(2)
%
% SEE ALSO:
%   `grid_structure`, `computeGeometry`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   opt = struct('maxArea',        10, ...
                'minAngle',       32, ...  % Degrees
                'segmentmarkers', [], ...
                'verbose',        false);
   opt = merge_options(opt, varargin{:});

   if isempty(opt.segmentmarkers),
      opt.segmentmarkers = int32(1 : size(edges, 1));
   end

   check_input(opt, points, edges);

   [tri, vor] = triangulate(opt, points, edges);

   % Construct output grid
   G = struct('cells'  , cells(tri, vor), ...
              'faces'  , faces(tri, vor), ...
              'nodes'  , nodes(tri),      ...
              'griddim', 2,               ...
              'type'   , { { mfilename } });
end

%--------------------------------------------------------------------------

function check_input(opt, points, edges)
   assert (numel(opt.segmentmarkers) == size(edges, 1), ...
           'There must be one segment identifier for each edge');

   assert (size(points, 2) == 2, 'Node coordinates must be planar (2D)');
   assert (size(edges , 2) == 2, 'Edges must connect exactly two nodes');

   assert (min(edges(:)) >  0, ...
           'Edges must refer to valid nodes (too small)');
   assert (max(edges(:)) <= size(points, 1), ...
           'Edges must refer to valid nodes (too large)');

   assert (opt.maxArea  >  0, 'Maximum area must be positive');
   assert (opt.minAngle >  0, 'Minumum angle must be positive');
   assert (opt.minAngle < 40, 'Minimum angle must be less than 40 degrees');
end

%--------------------------------------------------------------------------

function [tri, vor] = triangulate(opt, points, edges)
   % Domain boundary
   in = struct('pointlist'        , points .'          , ...
               'segmentlist'      , int32(edges .' - 1), ...
               'segmentmarkerlist', opt.segmentmarkers);

   options = ['Q2pzvAnei', ...
              sprintf('q%fa%f', opt.minAngle, opt.maxArea)];
   if opt.verbose,
      options = ['V', options(2:end)];
   end

   % Construct triangulation
   [tri, vor] = mex_triangle(in, options);
end

%--------------------------------------------------------------------------
% Computing cells.faces from triangle list requires some kind of sorting.

function c = cells(tri, vor)
   n  = size(tri.trianglelist, 2);
   cf = cellFaces(tri, vor);

   c  = struct('num'    , n,         ...
               'facePos', pos(3, n), ...
               'faces'  , cf(cf(:,1) ~= 0, 2));
end

%--------------------------------------------------------------------------

function f = faces(tri, vor)
   n = size(tri.edgelist, 2);

   f = struct('num'      , n,                                    ...
              'nodePos'  , pos(2, n),                            ...
              'nodes'    , double(reshape(tri.edgelist, [], 1)), ...
              'neighbors', double(vor.edgelist) .',              ...
              'tag'      , double(tri.edgemarkerlist) .');
end

%--------------------------------------------------------------------------

function n = nodes(tri)
   n = struct('num'   , size(tri.pointlist, 2), ...
              'coords', tri.pointlist .');
end

%--------------------------------------------------------------------------

function cf = cellFaces(tri, vor)
   e  = (1 : size(tri.edgelist, 2)) .';

   cf = double(vor.edgelist) .';
   cf = sortrows([ cf(:,1) , e ; ...
                   cf(:,2) , e ]);
end

%--------------------------------------------------------------------------

function v = pos(n, rpt)
   v = cumsum([1 ; repmat(n, [ rpt, 1 ])]);
end
