function G = initEclipseGrid(deck, varargin)
%Construct MRST grid from ECLIPSE GRID section.
%
% SYNOPSIS:
%   G = initEclipseGrid(deck)
%   G = initEclipseGrid(deck, 'pn1', pv1, ...)
%
% PARAMETERS:
%   deck - Raw input data in Deck form as defined by function
%          'readEclipseDeck'.
%
% OPTIONAL PARAMETERS:
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%            mapAxes - Apply the mapAxes transformation if present. 
%                      DEFAULT: false.
% 
%            In addition, the routine can pass on a number of paramters
%            to the function processGRDECL, which is used to construct
%            corner-point grids
%
% RETURNS:
%   G - Valid 'grid_structure'
%
% SEE ALSO:
%   `readEclipseDeck`, `grid_structure`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


   opt = struct('mapAxes', false);
   [opt,extra] = merge_options(opt, varargin{:});

   % -- Corner point grid -------------------------------------------------
   if all(isfield(deck.GRID, {'COORD', 'ZCORN'})) && ...
         isfield(deck.RUNSPEC, 'DIMENS')

      G = processGRDECL(deck.GRID, 'RepairZCORN', true, extra{:});

   % -- Tensor product grid -----------------------------------------------
   elseif all(isfield(deck.GRID, {'DXV', 'DYV', 'DZV'}))
      if isfield(deck.GRID, 'DEPTHZ'),
         depthz = deck.GRID.DEPTHZ;
      else
         depthz = zeros(prod(deck.RUNSPEC.DIMENS(1:2)+1), 1);
      end
      G = tensorGrid(cumsum([0; reshape(deck.GRID.DXV, [], 1)]), ...
                     cumsum([0; reshape(deck.GRID.DYV, [], 1)]), ...
                     cumsum([0; reshape(deck.GRID.DZV, [], 1)]), ...
                     'depthz', depthz);


   % -- Curvilinear logically Cartesian grid ------------------------------
   % -- FrontSim ----------------------------------------------------------
   elseif all(isfield(deck.GRID, {'COORDX', 'COORDY', 'COORDZ'})),
      G = buildCoordGrid(deck.GRID);


   % -- "Grid" ------------------------------------------------------------
   elseif all(isfield(deck.GRID, {'DX', 'DY', 'DZ'}))
      if isfield(deck.GRID, 'NNC')
         error('Non-neighboring connections not supported.');
      end

      nx = deck.RUNSPEC.DIMENS(1);
      ny = deck.RUNSPEC.DIMENS(2);
      nz = deck.RUNSPEC.DIMENS(3);

      dx = reshape(deck.GRID.DX, deck.RUNSPEC.DIMENS);
      dy = reshape(deck.GRID.DY, deck.RUNSPEC.DIMENS);
      dz = reshape(deck.GRID.DZ, deck.RUNSPEC.DIMENS);

      [dxv, n]= rlencode(dx,2);
      assert(n==ny, 'Only tensor-grid supported.');
      [dxv, n]= rlencode(dxv,3);
      assert(n==nz, 'Only tensor-grid supported.');

      [dyv, n]= rlencode(dy,1);
      assert(n==nx, 'Only tensor-grid supported.');
      [dyv, n]= rlencode(dyv,3);
      assert(n==nz, 'Only tensor-grid supported.');

      [dzv, n]= rlencode(dz,1);
      assert(n==nx, 'Only tensor-grid supported.');
      [dzv, n]= rlencode(dzv,2);
      assert(n==ny, 'Only tensor-grid supported.');


      if isfield(deck.GRID, 'TOPS'),
         if all(deck.GRID.TOPS(1:nx*ny) == deck.GRID.TOPS(1)),
            depthz = repmat(deck.GRID.TOPS(1), [(nx+1)*(ny+1), 1]);
         else
            error(['Block-centred grids are only supported ', ...
                   'for constant ''TOPS''.']);
         end
      else
         depthz = zeros([nx+1, ny+1]);
      end

      G = tensorGrid(cumsum([0; reshape(dxv, [], 1)]), ...
                     cumsum([0; reshape(dyv, [], 1)]), ...
                     cumsum([0; reshape(dzv, [], 1)]), ...
                     'depthz', depthz);
   else
      error('Grid not implemented')
   end

   if isfield(deck.GRID, 'MAPAXES') && opt.mapAxes
      for i = 1 : numel(G),
         G(i).nodes.coords(:,1:2) = ...
            mapAxes(G(i).nodes.coords(:,1:2), deck.GRID.MAPAXES);
      end
   end
end
