function src = evalBasisSource(G, weighting, rock)
%Evaluate synthetic basis function driving source term.
%
% SYNOPSIS:
%   src = evalBasisSource(G, w, rock)
%
% PARAMETERS:
%   G    - Grid data structure.
%   w    - Source weighting mode.  May be one of:
%            - Name (string) of particular weighting mode.  Must be one of
%                - 'perm' -- Weigh sources according to TRACE(K).
%                - 'poro'/'poros' --
%                            Weigh sources according to porosity.
%                - 'unit' -- Uniform sources.
%
%            - Array of already evaluated synthetic source terms.  This
%              array will be returned unchanged, and may be useful when
%              updating basis functions to account for compressible
%              effects.  The array is assumed to contain one non-negative
%              scalar value for each cell in the model.
%
%   rock - Rock data structure.  May be empty, but in the case of
%            - w == 'perm' -> must contain valid field 'rock.perm'.
%            - w == 'poro'/'poros'
%                          -> must contain valid field 'rock.poro'.
%
% RETURNS:
%   src  - Synthetic basis function source term for use with generator
%          functions 'evalBasisFunc' and 'evalWellBasis'.
%
% SEE ALSO:
%   `generateCoarseSystem`, `generateCoarseWellSystem`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


   if isnumeric(weighting) && numel(weighting) == G.cells.num,
      src = weighting;
   elseif ischar(weighting),
      ix = find(strcmp(weighting, {'perm', 'poro', 'poros', 'unit'}));

      if isempty(ix),
         error(msgid('BasisWeighting:Unknown'), ...
               'Basis function weighting strategy ''%s'' is unknown', ...
               weighting);
      end

      switch ix,
         case 1,
            if isempty(rock) || ~isfield(rock, 'perm'),
               error(msgid('RockData:Empty'), ...
                    'Rock data must contain valid field ''perm''.');
            end

            dim = size(G.nodes.coords, 2);
            K   = permTensor(rock, dim);

            assert (size(K,1) == G.cells.num, ...
                   ['Permeability must be defined in active cells',  ...
                    'only.\nGot %d tensors, expected %d (== number', ...
                    'of cells).'], size(K,1), G.cells.num);

            src = K(:, 1 : dim+1 : end) * ones([dim, 1]);  % == TRACE(K)

         case {2, 3},
            if isempty(rock) || ~isfield(rock, 'poro'),
               error(msgid('RockData:Empty'), ...
                     'Rock data must contain valid field ''poro''.');
            end

            src = rock.poro;

         case 4,

            src = ones([G.cells.num, 1]);

      end
   else
      error(msgid('BasisWeighting:NotSupported'), ...
           ['Basis function weighting does not conform to any ', ...
            'supported mode.']);
   end
end
