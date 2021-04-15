function op = setupOperatorsTPFA(G, rock, varargin)
% Set up helper structure for solvers based on automatic differentiation.
%
% SYNOPSIS:
%   op = setupOperatorsTPFA(G, rock)
%
% DESCRIPTION:
%   The automatic differentiation solvers rely on discrete operators for
%   divergence and gradient on the grid as well as a variety of derived
%   reservoir quantities such as transmissibility and pore volume. The
%   purpose of this function is to assemble all such quantities using a
%   standard two-point finite volume approxiation (TPFA).
%
% PARAMETERS:
%
%   G       - MRST grid given as a struct. See grid_structure.m for more
%             details.
%
%   rock    - Rock structure containing fields .perm and .poro with
%             approprioate dimensions for the grid. See makeRock for more
%             details.
%
% OPTIONAL PARAMETERS:
%   'trans'     - transmissibility for internal faces (if neighbors given)
%                 or for all faces (if neighbors are not given)
%
%   'neighbors' - neighbors for each internal face
%
%   'porv'      - pore volumes for all cells
%
% RETURNS:
%
%   op - Operators struct, with discrete operators and derived quantities:
%   
%      T_all - Transmissibilities for all interfaces, *including* (half)
%      transmissibilities for faces on the boundary. One value per
%      interface.
%
%      T - Transmissibilities for all internal interfaces. Internal
%      interfaces have a cell on both sides.
%
%      pv - Pore volumes. See function `poreVolume`. One value per cell.
%
%      faceAvg - (Function) For each interface, computes the average value
%      of a quantity defined in the cells. If a face is connecting two
%      cells, the faceAvg function will compute the arithmetic average of
%      the values in both cells.
%
%      M - Matrix used to compute face average
%     
%      internalConn - flag for internal connections (size G.faces.num)
%
%      Grad - (Function) Discrete gradient as function handle. Computes the
%      gradient on each interface via a first order finite difference
%      approximation using the values of the cells connected to the
%      face. Note that this discrete gradient does *not* divide by the
%      distance between the points.
%
%      C - Transfer matrix between cells and faces. Used to derive discrete
%      gradient and divergence operators.
%
%      Div - (Function) Discrete divergence. Integrates / sums up values on
%      the interfaces for all cells to give the (integrated) divergence per
%      cell.
%
%      AccDiv - (Function) adds accumulation term and discrete divergence term 
%
%      faceUpstr - (Function) Perform upstream weighting of values. Given a
%      set of cell wise values and a upstream flag for each interface, this
%      function will pick the values corresponding to the position in the
%      neighborship. I.e. if the flag is true for a given interface, it will
%      select the value in the FIRST cell connected to the interface
%      x(N(faceNo, 1)).  Otherwise, it will select the SECOND x(N(faceNo,
%      2)).  Typical usage is for upstream weighting of transported
%      quantities.
%
%      N - Neighborship structure. Will be number of interfaces by 2 in size
%      where N(ix, :) contains the cells connected to face number ix.
%
% SEE ALSO:
%   `computeTrans`, `processGRDECL`.

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

opt = struct('deck', [], 'neighbors', [], 'trans', [], 'porv', []);
opt = merge_options(opt, varargin{:});

if ~isempty(opt.deck)
   warning('DeckOption:Deprecated', ...
          ['The ''deck'' option is deprecated and inoperative. ', ...
           'It will be removed in a future version of MRST.']);
end

if isempty(opt.trans)
   [N, T, T_all, intInx] = grid_based_trans(G, rock, opt);
else
   [N, T, T_all, intInx] = user_provided_trans(G, opt.neighbors, opt.trans);
end

if has_grid_nncs(G)
   [N, T, T_all, intInx] = add_grid_nncs(N, T, T_all, intInx, G);
end

if any(T < 0)
   warn_negative_trans(T);
end

op = struct('T' , T,         'T_all',        T_all,  ...
            'N' , double(N), 'internalConn', intInx, ...
            'pv', pore_volume(G, rock, opt));

[op, nc, nf] = add_div_grad_operators(op, G, N);

op = add_face_average_operators(op, N, nf, nc);

% faceUpstr - as multiplication with matrix
op.faceUpstr = @(flag, x) faceUpstr(flag, x, N, [nf, nc]);

op.splitFaceCellValue = @(operators, flag, x) ...
   splitFaceCellValue(operators, flag, x, [nf, nc]);
end

%--------------------------------------------------------------------------

function [N, T, T_all, intInx] = grid_based_trans(G, rock, opt)
   assert(isempty(opt.neighbors))

   if isfield(G.faces, 'TRANS')
      T = G.faces.TRANS;
   else
      % half-trans -> trans and reduce to interior
      T = getFaceTransmissibility(G, rock);
   end

   N      = G.faces.neighbors;
   intInx = all(N ~= 0, 2);

   T_all = T;
   N     = N(intInx, :);
   T     = T(intInx);
end

%_-------------------------------------------------------------------------

function [N, T, T_all, intInx] = user_provided_trans(G, N, T)
   % transmissibility is given as input
   if isempty(N)
      [N, intInx, n_if] = internal_grid_neighbours(G);
   else
      [N, intInx, n_if] = explicit_cell_pairs(N, G);
   end

   if numel(T) == n_if
      % Internal interface transmissibility
      T_all = zeros(size(intInx));
      T_all(intInx) = T;

   else
      % All transmissibilities given
      assert(numel(T) == numel(intInx));

      T_all = T;
      T     = T(intInx);
   end
end

%--------------------------------------------------------------------------

function [N, T, T_all, intInx] = add_grid_nncs(N, T, T_all, intInx, G)
   cells = G.nnc.cells;
   trans = G.nnc.trans;

   i = ~ (any(cells == 0, 2) | ...
          ismember(sort(cells, 2), sort(N, 2), 'rows'));

   N      = [ N      ; cells(i, :) ];
   T      = [ T      ; trans(i) ];
   T_all  = [ T_all  ; trans(i) ];
   intInx = [ intInx ; true([sum(i), 1]) ];
end

%--------------------------------------------------------------------------

function pv = pore_volume(G, rock, opt)
   pv = opt.porv;

   if isempty(pv)
      if isfield(G.cells, 'PORV')
         pv = G.cells.PORV;
      else
         pv = poreVolume(G, rock);
      end

      numzeropv = sum(~ (pv > 0));
      if numzeropv > 0
         warning('ZeroPV:InActiveCells', ...
                ['I computed zero pore volumes in %d cells.  ', ...
                 'Consider adjusting ''poro''/''ntg'' fields ', ...
                 'or the grid itself.'], numzeropv);
      end
   end
end

%--------------------------------------------------------------------------

function [op, nc, nf] = add_div_grad_operators(op, G, N)
   [nc, nf] = deal(numel(op.pv), size(N, 1));

   assert(nc == G.cells.num, ...
          'Dimension mismatch between grid and supplied pore-volumes.');

   % C - (transpose) divergence matrix
   C  = sparse([(1:nf)'; (1:nf)'], N, ones(nf, 1) * [1, -1], nf, nc);

   op.C    = C;
   op.Grad = @(x) -C*x;

   if nf == 0
      % We have zero faces. Account for Matlab's preference for
      % reducing expressions of type a + [] to [].
      op.AccDiv = @(acc, flux) acc;
      op.Div = @(x) zeros(nc, 1);
   else
      op.AccDiv = @(acc, flux) acc + C'*flux;
      op.Div = @(x) C'*x;
   end
end

%--------------------------------------------------------------------------

function op = add_face_average_operators(op, N, nf, nc)
   M = sparse((1 : nf) .' * [1, 1], N, repmat(0.5, [nf, 2]), nf, nc);

   op.M = M;
   op.faceAvg = @(x) M*x;
end

%--------------------------------------------------------------------------

function warn_negative_trans(T)
   nneg = sum(T < 0);

   warning('Transmissibility:Negative', ...
          ['Negative transmissibilities detected for %d/%d (%.2e %%) ', ...
           'connections.'], nneg, numel(T), 100 * nneg / numel(T));
end

%--------------------------------------------------------------------------

function tf = has_grid_nncs(G)
   tf = isfield(G, 'nnc') && isstruct(G.nnc) && ...
      all(isfield(G.nnc, {'cells', 'trans'}));
end

%--------------------------------------------------------------------------

function [N, intInx, n_if] = internal_grid_neighbours(G)
   % Get neighbors for internal faces from grid.
   N      = double(G.faces.neighbors);
   intInx = all(N ~= 0, 2);

   N    = N(intInx, :);
   n_if = sum(intInx);
end

%--------------------------------------------------------------------------

function [N, intInx, n_if] = explicit_cell_pairs(N, G)
   % neighbors are given
   intInx = all(N ~= 0, 2);
   n_if   = sum(intInx);

   if isfield(G, 'faces')
      % Try to match given interfaces to actual grid.
      intInxGrid = all(G.faces.neighbors ~= 0, 2);

      if sum(intInxGrid) == n_if
         % Given neighbors correspond to internal interfaces
         intInx = intInxGrid;
      end
   end
end
