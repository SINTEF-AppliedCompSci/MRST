function m = computeTranMult(G, rockprop)
%Compute transmissibility multipliers.
%
% SYNOPSIS:
%   mult = computeTranMult(G, rockprop)
%
% PARAMETERS:
%   G    - Valid grid structure.  Must be a logically Cartesian grid such
%          as the grids created by means of functions 'processGRDECL' or
%          'cartGrid'.  Moreover, the half-face tags (G.cells.faces(:,2))
%          must not contain values outside the set (1 : 6).
%
%   rockprop -
%          Rock property data structure.  May contain ECLIPSE multiplier
%          keywords such as 'MULTX', 'MULTY', or 'MULTZ-' (matched case
%          insensitively).  Often corresponds to the 'deck.GRID' data
%          structure of function 'readEclipseDeck'.
%
% RETURNS:
%   mult - Half-face transmissibility multipliers.  An empty array (i.e.,
%          mult = []) if the rock property data structure contains no
%          multiplier keyword fields.  Otherwise (i.e., when there are
%          multiplier keywords present), a SIZE(G.cells.faces,1)-by-1
%          DOUBLE array such that 'mult(i)' is a transmissibility
%          multiplier for half-face (half-contact) 'i'.
%
% NOTE:
%   A multiplier keyword, if present, must contain one datum for each cell
%   in the global grid (PROD(G.cartDims) cells).  Values in inactive cells
%   (e.g., due to explicit deactivation or pinch/collapse) are ignored.
%
%   We underscore the fact that these multiplier values are designed to go
%   with the half-face (one-sided) transmissibility values derived from
%   function 'computeTrans'.  In particular, the multiplicative factors
%   must be included BEFORE computing the harmonic average of the one-sided
%   transmissibilities.
%
% EXAMPLE:
%   %% Compute connection transmissibilities on an ECLIPSE-type model
%   %
%   % 1) Compute background, one-sided transmissibilities
%   htrans = computeTrans(G, rock);
%
%   % 2) Compute multipliers for one-sided transmissibilities
%   tmult = computeTranMult(G, deck.GRID);
%
%   if ~ isempty(tmult),
%      % 2b) Include multipliers for one-sided transmissibilities
%      htrans = htrans .* tmult;
%   end
%
%   % 3) Compute final background (static) transmissibilities
%   trans = 1 ./ accumarray(G.cells.faces(:,1), ...
%                           1 ./ htrans,        ...
%                           [G.faces.num, 1]);
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

   if isempty(rockprop)
       % No input given
       m = [];
       return
   end

   if size(G.cells.faces, 2) < 2
      error(['Function ''%s'' is only supported in corner-point ', ...
             'grids which identify cardinal directions.'], mfilename);
   end

   multkw     = fieldnames(rockprop);
   kw_pattern = '^mult([xyz][-_]?)$';

   i = ~cellfun(@isempty, regexpi(multkw, kw_pattern, 'once'));

   if ~any(i)
      m = [];
   else
      multkw = multkw(i);
      dir_id = identify_directions();

      ix   = [G.cells.faces(:,1), ...
              rldecode(G.cells.indexMap, diff(G.cells.facePos))];

      subs = [];
      vals = [];

      n_exptd = prod(G.cartDims);

      j = false([max(G.cells.faces(:,2)), 1]);
      for kw = reshape(multkw, 1, [])
         prop = reshape(rockprop.(kw{1}), [], 1);

         if numel(prop) ~= n_exptd
            assert (numel(prop) < n_exptd, ...
                   ['Multiplier keyword ''%s'' specifies more values ', ...
                    '(%d) than there are cells in global grid (%d)\n'], ...
                    regexprep(kw{1}, '_', '-'), numel(prop), n_exptd);

            prop = [prop; ones([n_exptd - numel(prop), 1])];  %#ok
         end

         prop(~isfinite(prop)) = 1; % Suitable default value for multiplier

         d    = dir_id(regexprep(kw{1}, kw_pattern, '$1', 'ignorecase'));

         j(d) = true;
         i    = ix(j(G.cells.faces(:,2)), :);
         j(d) = false;

         subs = [subs;              i(:, 1)         ];     %#ok
         vals = [vals; reshape(prop(i(:, 2)), [], 1)];     %#ok
      end

      m = accumarray(subs, vals, [G.faces.num, 1], @prod, 1);
      m = m(ix(:, 1));
   end
end

%--------------------------------------------------------------------------

function id = identify_directions()
   s  = struct('x_', 1, 'i_', 1, 'x', 2, 'i', 2, ...
               'y_', 3, 'j_', 3, 'y', 4, 'j', 4, ...
               'z_', 5, 'k_', 5, 'z', 6, 'k', 6);
   id = @(d) s.(lower(regexprep(d, '\W', '_')));
end
