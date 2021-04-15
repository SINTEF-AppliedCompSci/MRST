function G = mprocessGRDECL(grdecl, varargin)
%Compute grid topology and geometry from pillar grid description.
%
% SYNOPSIS:
%   G = mprocessGRDECL(grdecl)
%   G = mprocessGRDECL(grdecl, 'pn1', pv1, ...)
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined by function
%            'readGRDECL', with fields COORDS, ZCORN and, possibly, ACTNUM.
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              Verbose -- Whether or not to display progress information
%                         Logical.  Default value: Verbose = false.
%
%              Tolerance --
%                         Minimum distinguishing vertical distance for
%                         points along a pillar.  Specifically, two points
%                         (x1,y1,z1) and (x2,y2,z2) are considered separate
%                         only if ABS(z2 - z1) > Tolerance.
%                         Non-negative scalar.
%                         Default value: Tolerance = 0.0 (distinguish all
%                         points along a pillar whose z coordinate differ
%                         even slightly).
%
%              CheckGrid --
%                         Whether or not to perform basic consistency
%                         checks on the resulting grid.
%                         Logical.  Default value: CheckGrid = true.
%
%              SplitDisconnected --
%                         Whether or not to split disconnected grid
%                         components into separate grids/reservoirs.
%                         Logical.  Default value: SplitDisconnected=true.
%
% RETURNS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
% EXAMPLE:
%   G = mprocessGRDECL(readGRDECL('small.grdecl'));
%   plotGrid(G); view(10,45);
%
% SEE ALSO:
%   `readGRDECL`, `processGRDECL`

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

   opt = struct('Verbose', mrstVerbose, 'Tolerance', 0.0, ...
                'CheckGrid', true, 'SplitDisconnected', true, ...
                'PreserveCpNodes', false);
   opt = merge_options(opt, varargin{:});

   assert (~opt.PreserveCpNodes, ...
          ['Function %s does not currently implement the ', ...
           '''PreserveCpNodes'' option'], mfilename);

   if isfield(grdecl, 'ACTNUM') && ...
         ~isa(grdecl.ACTNUM, 'int32')

      grdecl.ACTNUM = int32(grdecl.ACTNUM);
   end

   G = processgrid_mex(grdecl,opt.Tolerance);

   % Guarantee type DOUBLE for indices.  Needed, for instance, if the grid
   % is used in constructing linear operators in the AD framework (which
   % pass such index arrays to function SPARSE).
   G = double_indices(G);

   if opt.CheckGrid
      assert (all(diff(G.cells.facePos) > 3), ...
              'All cells must have at least four faces');

      assert (all(diff(G.faces.nodePos) > 2), ...
              'All faces must have at least three nodes');

      assert (all(all(isfinite(G.nodes.coords))), ...
              'All nodes must have all finite coordinates');
   end

   if opt.SplitDisconnected
      G = splitDisconnectedGrid(G, 'Verbose', false);
   end

   [ G.type    ] = deal({ mfilename });
   [ G.griddim ] = deal(3);

%{
   if isfield(grdecl, 'MAPAXES'),
      for i = 1 : numel(G),
         G(i).nodes.coords(:,1:2) = ...
            mapAxes(G(i).nodes.coords(:,1:2), grdecl.MAPAXES);
      end
   end
%}
end

%--------------------------------------------------------------------------

function G = double_indices(G)
   G.faces.neighbors = double(G.faces.neighbors);
   G.faces.nodes     = double(G.faces.nodes);
   G.faces.nodePos   = double(G.faces.nodePos);

   G.cells.faces     = double(G.cells.faces);
   G.cells.facePos   = double(G.cells.facePos);

   if isfield(G.cells, 'indexMap')
      G.cells.indexMap = double(G.cells.indexMap);
   end
end
