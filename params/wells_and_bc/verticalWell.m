function W = verticalWell(W, G, rock, I, varargin)
%Insert a vertical well into the simulation model.
%
% SYNOPSIS:
%   Either 1)
%      W = verticalWell(W, G, rock, I, J, K)
%      W = verticalWell(W, G, rock, I, J, K, 'pn1', pv1, ...)
%
%   Or 2)
%      W = verticalWell(W, G, rock, I,    K)
%      W = verticalWell(W, G, rock, I,    K, 'pn1', pv1, ...)
%
% PARAMETERS:
%   W    - Input well structure defining existing wells.  Pass an empty
%          structure if there are no previously defined wells.  This
%          structure is updated on output.
%
%   G    - Grid data structure.  Must contain valid field
%          `G.cells.centroids`.  Call function `computeGeometry` to obtain
%          these values.
%
%   rock - Rock data structure.  Must contain valid field `rock.perm`.
%
%   I,J  - Horizontal location of well heel. Two possible modes,
%
%          Mode 1:
%            `I` and `J` are Cartesian indices.  Specifically,
%            `I` is the index along the first logical direction while
%            `J` is the index along the second logical direction.
%
%            This mode is only supported in grids which have an underlying
%            Cartesian (logical) structure such as purely Cartesian grids or
%            corner-point grids.
%
%          Mode 2:
%            `I` is the cell index `(1 <= I <= G.cells.num)` of the
%            *top-most* cell in the column through which the vertical well
%            will be completed.  `J` must not be present.
%
%            This mode is supported for logically Cartesian grids containing
%            a three-component field `G.cartDims` or for otherwise layered
%            grids which contain the fields `G.numLayers` and `G.layerSize`.
%
%   K    - A vector of layers in which this well should be completed.
%          An empty layer vector (i.e., `isempty(K)`), is replaced by::
%
%               K = 1 : num_layers
%
%          In other words, `isempty(K)` implies completion in ALL layers.
%
% OPTIONAL PARAMETERS:
%   'Any' - All options supported by function `addWell` are supported in
%           `verticalWell`.  Any direction specifications (i.e., option
%           'Dir') are ignored and replaced by 'z'.
%
% RETURNS:
%   W - Updated well structure.
%
% EXAMPLE:
%   G = cartGrid([60, 220, 85], [60, 220, 85].*[20, 10, 2].*ft);
%   G = computeGeometry(G);
%   rock.perm = repmat([500, 500, 50].*milli*darcy(), [G.cells.num, 1]);
%   W = struct([]);
%   W = verticalWell(W, G, rock,  1,   1, (1:85),     ...
%                    'Type', 'bhp', 'Val', 300*barsa, ...
%                    'Radius', 0.125, 'Name', 'P1');
%   W = verticalWell(W, G, rock, 60,   1, (1:85),     ...
%                    'Type', 'bhp', 'Val', 300*barsa, ...
%                    'Radius', 0.125, 'Name', 'P2');
%   W = verticalWell(W, G, rock,  1, 220, (1:85),     ...
%                    'Type', 'bhp', 'Val', 300*barsa, ...
%                    'Radius', 0.125, 'Name', 'P3');
%   W = verticalWell(W, G, rock, 60, 220, (1:85),     ...
%                    'Type', 'bhp', 'Val', 300*barsa, ...
%                    'Radius', 0.125, 'Name', 'P4');
%   W = verticalWell(W, G, rock, 30, 110, (1:85),     ...
%                    'Type', 'bhp', 'Val', 500*barsa, ...
%                    'Radius', 0.125, 'Name', 'I1',   ...
%                    'Comp_i', [0.5, 0.5, 0]);
%
% SEE ALSO:
%   `grid_structure`, `makeLayeredGrid`, `addWell`.

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


assert (isnumeric(I) && numel(I) == 1);

%--------------------------------------------------------------------------
% Determine call mode. ----------------------------------------------------
%
if mod(nargin, 2) == 0
   % W = verticalWell(W, G, rock, I, J, K, ...)
   %
   assert (isfield(G, 'cartDims'));
   layerSize = prod(G.cartDims(1:2));
   if G.griddim == 3
      numLayers = G.cartDims(3);
   else
      numLayers = 1;
   end

   J = varargin{1}; varargin = varargin(2 : end);
   assert (isnumeric(J) && numel(J) == 1);

   assert ((1 <= I) && (I <= G.cartDims(1)), 'I-position out of bounds');
   assert ((1 <= J) && (J <= G.cartDims(2)), 'J-position out of bounds');

   I = I + (J - 1)*G.cartDims(1);    clear J
else
   % W = verticalWell(W, G, rock, I, K, ...)
   %
   if isfield(G, 'layerSize') && isfield(G, 'numLayers')
      % G presumably a 'makeLayeredGrid'.
      layerSize = G.layerSize;
      numLayers = G.numLayers;
   elseif isfield(G, 'cartDims')
      % Logically Cartesian grid.
      layerSize = prod(G.cartDims(1:2));
      if G.griddim == 3
         numLayers = G.cartDims(3);
      else
         numLayers = 1;
      end
   else
      error(msgid('GridType:NotLayered'), ...
            'Grid type does not appear to support a layered structure');
   end
end

%--------------------------------------------------------------------------
% Extract layer vector (K). -----------------------------------------------
%
K = varargin{1}; varargin = varargin(2 : end);
if isempty(K)
   % Empty 'K'.  Complete in all layers.
   K = (1 : numLayers) .';
end
assert ((1 <= min(K)) && (max(K) <= numLayers));

%--------------------------------------------------------------------------
% Extract active completions. ---------------------------------------------
%
act                   = zeros([numLayers * layerSize, 1]);
act(G.cells.indexMap) = 1 : numel(G.cells.indexMap);

wc = act(I + (K - 1).*layerSize);
i  = wc > 0;

%--------------------------------------------------------------------------
% Define resulting well. --------------------------------------------------
%
W = addWell(W, G, rock, wc(i), varargin{:}, 'dir', 'z');
