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

   opt = struct('mapAxes',      false, ...
                'removeZeroPV', false, ...
                'useMex',       false);
   [opt, extra] = merge_options(opt, varargin{:});
   % -- Corner point grid -------------------------------------------------
   if all(isfield(deck.GRID, {'COORD', 'ZCORN'})) && ...
         isfield(deck.RUNSPEC, 'DIMENS')
     if opt.useMex
         % Use MEX accelerated version
         G = mprocessGRDECL(deck.GRID, extra{:});
     else
         % Use Matlab version
         G = processGRDECL(deck.GRID, 'RepairZCORN', true, extra{:});
     end

   % -- Tensor product grid -----------------------------------------------
   elseif all(isfield(deck.GRID, {'DXV', 'DYV', 'DZV'}))
      if isfield(deck.GRID, 'DEPTHZ')
         depthz = deck.GRID.DEPTHZ;
      else
         depthz = zeros(prod(deck.RUNSPEC.DIMENS(1:2)+1), 1);
      end
      G = tensorGrid(cumsum([0; reshape(deck.GRID.DXV, [], 1)]), ...
                     cumsum([0; reshape(deck.GRID.DYV, [], 1)]), ...
                     cumsum([0; reshape(deck.GRID.DZV, [], 1)]), ...
                     'depthz', depthz);

      if isfield(deck.GRID, 'ACTNUM')
         G = extractSubgrid(G, deck.GRID.ACTNUM > 0);
      end

   % -- Curvilinear logically Cartesian grid ------------------------------
   % -- FrontSim ----------------------------------------------------------
   elseif all(isfield(deck.GRID, {'COORDX', 'COORDY', 'COORDZ'}))
      if exist('buildCoordGrid', 'file')
         % This version of MRST has experimental support for COORD{X,Y,X}.
         G = buildCoordGrid(deck.GRID);
      else
         error('GridType:Unsupported', ...
               'MRST does not support the COORD{X,Y,Z} grid type');
      end

   % -- "MRST-output NNC grid" --------------------------------------------
   elseif is_nnc_1d_grid(deck)
      nc = max(deck.GRID.cartDims);
      e = deck.EDIT;
      nnc = deck.GRID.NNC;
      N = nnc(:, [1, 4]);
      T = nnc(:, 7);
      G = cartGrid([nc, 1, 1]);
      G.cells.faces = [1; 1];
      G.cells.eMap = ':';
      G.faces.num = 0;
      rock = makeRock(G, deck.GRID.PERMX, deck.GRID.PORO);
      G = simGridTPFA(G, rock, 'neighbors', N, 'porv', e.PORV, 'depth', e.DEPTH, 'trans', T);

   % -- "Grid" ------------------------------------------------------------
   elseif is_delta_grid(deck)
      nx = deck.RUNSPEC.DIMENS(1);
      ny = deck.RUNSPEC.DIMENS(2);
      nz = deck.RUNSPEC.DIMENS(3);

      dx = getDeltas(deck, 1);
      dy = getDeltas(deck, 2);
      dz = getDeltas(deck, 3);

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


      if isfield(deck.GRID, 'TOPS')
         if all(deck.GRID.TOPS(1:nx*ny) == deck.GRID.TOPS(1))
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

      if isfield(deck.GRID, 'ACTNUM')
         G = extractSubgrid(G, deck.GRID.ACTNUM > 0);
      end

   else
      error('Grid not implemented')
   end

   if isfield(deck.GRID, 'NNC')
      G = merge_nnc(G, processNNC(G, deck.GRID.NNC));
   end

   if isfield(deck.GRID, 'MAPAXES') && opt.mapAxes
      for i = 1:numel(G)
         G(i).nodes.coords(:,1:2) = ...
            mapAxes(G(i).nodes.coords(:,1:2), deck.GRID.MAPAXES);
      end
   end

   if opt.removeZeroPV
       mask = true(prod(deck.GRID.cartDims), 1);
       if isfield(deck.GRID, 'PORO')
           mask = mask & deck.GRID.PORO > 0;
       end
       if isfield(deck.GRID, 'PORV')
           mask = mask & deck.GRID.PORV > 0;
       end
       for i = 1:numel(G)
           G(i) = extractSubgrid(G(i), mask(G(i).cells.indexMap));
       end
   end
end

function is_grid = is_delta_grid(deck)
    is_grid = all(isfield(deck.GRID, {'DX', 'DY', 'DZ'}) | ...
                  isfield(deck.GRID, {'DXV', 'DYV', 'DZV'}));
end

function is_nnc = is_nnc_1d_grid(deck)
    is_nnc = all(isfield(deck.GRID, {'DX', 'DY', 'DZ'})) && ...
        isfield(deck.GRID, 'cartDims') && ...
        prod(deck.GRID.cartDims) == max(deck.GRID.cartDims) && ...
        isfield(deck.GRID, 'NNC');
end

function dc = getDeltas(deck, index)
xyz = 'XYZ';
duv = ['D', xyz(index), 'V'];
du = ['D', xyz(index)];
dims = deck.RUNSPEC.DIMENS;

if isfield(deck.GRID, duv)
    % Single vector, expand to nx by ny by nz format
    arg1 = {dims(1), dims(2), dims(3)};
    arg2 = {1, 1, 1};
    arg1{index} = 1;
    arg2{index} = [];
    dxv = reshape(deck.GRID.(duv), arg2{:});
    dc = repmat(dxv, arg1{:});
else
    % We were given deltas directly
    dc = reshape(deck.GRID.(du), dims);
end
end

%--------------------------------------------------------------------------

function G = merge_nnc(G, nnc)
   if ~isfield(G, 'nnc')
      G.nnc = nnc;
   else
      G.nnc = do_merge_nnc(G.nnc, nnc);
   end
end

%--------------------------------------------------------------------------

function nnc = do_merge_nnc(nnc, new)
   i = ~ismember(sort(new.cells, 2), sort(nnc.cells, 2), 'rows');

   if any(i)
      nnc.cells = [ nnc.cells ; new.cells(i,:) ];
      nnc.trans = [ nnc.trans ; new.trans(i) ];
   end
end
