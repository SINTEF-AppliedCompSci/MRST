function [K, r, c] = permTensor(rock, dim)
%Expand permeability tensor to full format.
%
% SYNOPSIS:
%    K        = permTensor(rock, dim)
%   [K, r, c] = permTensor(rock, dim)
%
% PARAMETERS:
%   rock - Rock data structure with valid field `perm`--one value of the
%          permeability tensor for each cell (nc) in the discretised model.
%
%          The field `rock.perm` may have ONE column for a scalar
%          (isotropic) permeability in each cell, TWO or THREE columns for
%          a diagonal permeability in each cell (in two or three space
%          dimensions, respectively) and THREE or SIX columns for a
%          symmetric, full tensor permeability.  In the latter case, each
%          cell gets the permeability tensor::
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
%   dim  - Number of space dimensions in the discretised model.  Must be
%          two or three.
%
% RETURNS:
%   K    - Expanded permeability tensor.  An `nc` by `dim^2` array of
%          permeability values as described above.
%
%   r, c - Row- and column indices from which the two-form
%
%                (x, Ky) = \sum_r \sum_c x(r) * K(r,c) * y(c)
%
%          may be easily evaluated by a single call to SUM (see EXAMPLE).
%          OPTIONAL.  Only returned if specifically requested.
%
% EXAMPLE:
%   % Compute n'Kg gravity term on each cell face (half contact).
%   %
%   cellNo = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2).';
%   sgn    = 2*double(G.faces.neighbors(G.cells.faces(:,1), 1) == cellNo) - 1;
%   % cfn  = cell-face normal = *outward* face normals on each cell.
%   cfn    = bsxfun(@times, G.faces.normals(G.cells.faces(:,1), :), sgn);
%
%   dim    = size(G.nodes.coords, 2);
%   [K, r, c] = permTensor(rock, dim);
%
%   g   = gravity();   g = reshape(g(1 : dim), 1, []);
%   nKg = sum(cfn(:,r) .* bsxfun(@times, K(cellNo,:), g(c)), 2);
%
% SEE ALSO:
%   `computeTrans`.

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

   if isempty(rock) || ~isfield(rock, 'perm')
      error(msgid('Rock:Empty'), ...
            '''rock'' structure must contain valid field ''perm''.')
   end

   [nc, nk] = size(rock.perm);
   if dim == 1
      assert (nk == 1, ...
             ['A one dimensional domain does not support ', ...
              'multi-component tensors.']);

      K = rock.perm;
      r = 1;
      c = 1;

   elseif dim == 2
      K = symmetric_2D(rock, nk, nc);

      K = K(:, [1, 2, 2, 3]);
      r =      [1, 1, 2, 2] ;
      c =      [1, 2, 1, 2] ;

   elseif dim == 3
      K = symmetric_3D(rock, nk, nc);

      K = K(:, [1, 2, 3, 2, 4, 5, 3, 5, 6]);
      r =      [1, 1, 1, 2, 2, 2, 3, 3, 3] ;
      c =      [1, 2, 3, 1, 2, 3, 1, 2, 3] ;

   else
      error(msgid('PermDim:Unsupported'), ...
            '%d space dimensions are unsupported at this time.', dim);
   end
end

%--------------------------------------------------------------------------

function K = symmetric_2D(rock, nk, nc)
   switch nk
      case 1
         % Isotropic.

         K = rock.perm * [1, 0, 1];

      case 2
         % Diagonal.

         K = [rock.perm(:,1), zeros([nc,1]), rock.perm(:,2)];

      case 3
         % Full, symmetric, 2-by-2.

         K = rock.perm;

      otherwise
         error(msgid('PermDim:Unsupported'), ...
              ['%d-component permeability value is not supported ', ...
               'in two space dimensions.'], nk);
   end
end

%--------------------------------------------------------------------------

function K = symmetric_3D(rock, nk, nc)
   switch nk
      case 1
         % Isotropic.

         K = rock.perm * [1, 0, 0, 1, 0, 1];

      case 3
         % Diagonal.

         K = [rock.perm(:,1), zeros([nc,2]), ...
              rock.perm(:,2), zeros([nc,1]), ...
              rock.perm(:,3)];

      case 6
         % Full, symmetric, 3-by-3.

         K = rock.perm;

      otherwise
         error(msgid('PermDim:Unsupported'), ...
              ['%d-component permeability value is not supported ', ...
               'in three space dimensions.'], nk);
   end
end
