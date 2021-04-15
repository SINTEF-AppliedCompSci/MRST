function rock = grdecl2Rock(grdecl, varargin)
%Extract rock properties from input deck (grdecl structure)
%
% SYNOPSIS:
%   rock = grdecl2Rock(grdecl)
%   rock = grdecl2Rock(grdecl, actmap)
%
% PARAMETERS:
%   grdecl - Raw input data in `GRDECL` form (typically corresponding to
%            values from the `GRID` section of an ECLIPSE or FrontSim input
%            deck.
%
%   actmap - Active-to-global cell mapping.  An n-by-1 array of global cell
%            indices.  Specifically, actmap(i) is the linearized Cartesian
%            index of the global cell (from the nx-by-ny-by-nz bounding
%            box) corresponding to active cell `i`.
%
%            OPTIONAL.  Default value corresponds to treating all input
%            values as corresponding to active grid cells.
%
% RETURNS:
%   rock - Rock data stucture suitable for passing to function
%          `permTensor`.  Specifically, function `grdecl2Rock` assigns the
%          permeability tensor values, `rock.perm`.  Moreover, the
%          following data values will be assigned if present in the input
%          deck:
%            - poro -- Porosity values.  Corresponds to `PORO` keyword.
%            - ntg  -- Net-to-gross values.  Corresponds to `NTG` keyword.
%
%          If, additionally, the input 'grdecl' structure contains
%          transmissibility multiplier data (keywords 'MULTX', 'MULTY',
%          'MULTZ', 'MULTX-', 'MULTY-', 'MULTZ-'), then those values will
%          be stored as individual fields in a substructure 'multipliers'.
%          The field names are lower case without 'MULT' and hypens are
%          replaced by underscores.  For instance, the 'MULTX-' data will
%          be stored as
%
%             - rock.multipliers.x_
%
%          Furthermore, if the input 'grdecl', structure contains fault
%          multiplier data (both of the keywords 'FAULTS' and 'MULTFLT'),
%          then those values will be copied verbatim to a substructure
%          'faultdata' and retain their original field names, converted to
%          lower case, in this substructure.
%
% NOTE:
%   Function `grdecl2Rock` only extracts the raw data values from the input
%   vectors.  The caller will have to perform any required unit conversion
%   separately, possibly aided by functions `convertInputUnits` and
%   `convertFrom`.
%
%   For backwards compatibility, the caller may pass a valid grid_structure
%   as the second parameter (`actmap`) of function `grdecl2Rock`.  In this
%   case, the field `grid_structure.cells.indexMap` (if present and
%   non-empty) will be used in place of the `actmap` array described above.
%
% SEE ALSO:
%   `grid_structure`, `permTensor`, `convertInputUnits`, `convertFrom`.

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

   % Check consistency.
   if ~consistent(grdecl)
      error(msgid('Tensor:Inconsistent'), ...
            'Input tensor is not structurally consistent.');
   end

   % Permeability tensor.
   [rock.perm, actmap] = getTensor(grdecl, varargin{:});
   rock.perm = rock.perm(actmap, :);

   extract = @(fld) reshape(grdecl.(fld)(actmap), [], 1);

   % Other properties.
   rockprop = {'PORO', 'NTG', 'ROCKTYPE'};
   for fld = reshape(rockprop(isfield(grdecl, rockprop)), 1, [])
      rock.(lower(fld{1})) = extract(fld{1});
   end

   cname = @(kw) lower(kw{1}(5 : end));
   mult  = strcat('MULT', {'X', 'Y', 'Z'});
   mult  = [ mult, strcat(mult, '_') ];
   for prop = reshape(mult(isfield(grdecl, mult)), 1, [])
      rock.multipliers.(cname(prop)) = extract(prop{1});
   end

   if all(isfield(grdecl, {'FAULTS', 'MULTFLT'}))
      rock.faultdata = struct('faults' , { grdecl.FAULTS }, ...
                              'multflt', { grdecl.MULTFLT });
   end
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function b = consistent(grdecl)
   k = false([3, 3]);

   k(1,1) = isfield(grdecl, 'PERMX' ) || isfield(grdecl, 'PERMXX') || isfield(grdecl, 'PERMH');
   k(1,2) = isfield(grdecl, 'PERMXY') || isfield(grdecl, 'PERMYX');
   k(1,3) = isfield(grdecl, 'PERMXZ') || isfield(grdecl, 'PERMZX');

   k(2,1) = isfield(grdecl, 'PERMYX') || isfield(grdecl, 'PERMXY');
   k(2,2) = isfield(grdecl, 'PERMY' ) || isfield(grdecl, 'PERMYY') || isfield(grdecl, 'PERMH');
   k(2,3) = isfield(grdecl, 'PERMYZ') || isfield(grdecl, 'PERMZY');

   k(3,1) = isfield(grdecl, 'PERMZX') || isfield(grdecl, 'PERMXZ');
   k(3,2) = isfield(grdecl, 'PERMZY') || isfield(grdecl, 'PERMYZ');
   k(3,3) = isfield(grdecl, 'PERMZ' ) || isfield(grdecl, 'PERMZZ');

   b = any(k(:));
   b = b && (k(1,1) || ~any([k(1,2), k(2,1), k(1,3), k(3,1)]));
   b = b && (k(2,2) || ~any([k(2,1), k(1,2), k(2,3), k(3,2)]));
   b = b && (k(3,3) || ~any([k(3,1), k(1,3), k(3,2), k(2,3)]));
end

%--------------------------------------------------------------------------

function [perm, actmap] = getTensor(grdecl, varargin)
   [vals, comp] = tensorValues(grdecl);

   actmap = (1 : size(vals,1)) .';
   if nargin > 1
      if isnumeric(varargin{1})
         actmap = varargin{1};
      elseif isstruct(varargin{1})          && ...
            isfield  (varargin{1}, 'cells') && ...
            isstruct (varargin{1}.  cells ) && ...
            isfield  (varargin{1}.  cells , 'indexMap')
         % varargin{1} is likely to be a grid_structure.
         % Treat as such and extract its 'indexMap' field.
         actmap = varargin{1}.cells.indexMap;
      end
   end

   [i, j] = find(comp > 1);
   if all(i == j)
      % Diagonal.
      assert (numel(i) == 3);
   else
      % Full, symmetric.  Use upper triangular part.
      assert (numel(i) == 9);
      i = [1, 1, 1, 2, 2, 3] .';
      j = [1, 2, 3, 2, 3, 3] .';
   end

   comp = comp(sub2ind([3, 3], i, j));
   if numel(comp) == 3 && all(diff(comp) == 0)
      % Return only single-component tensor when isotropic.
      comp = comp(1);
   end

   perm = vals(:, comp);
end

%--------------------------------------------------------------------------

function [vals, comp] = tensorValues(grdecl)
   nc = 0;
   if nc == 0 && isfield(grdecl, 'PERMX' ), nc = numel(grdecl.PERMX ); end
   if nc == 0 && isfield(grdecl, 'PERMY' ), nc = numel(grdecl.PERMY ); end
   if nc == 0 && isfield(grdecl, 'PERMH' ), nc = numel(grdecl.PERMH ); end
   if nc == 0 && isfield(grdecl, 'PERMZ' ), nc = numel(grdecl.PERMZ ); end
   if nc == 0 && isfield(grdecl, 'PERMXX'), nc = numel(grdecl.PERMXX); end
   if nc == 0 && isfield(grdecl, 'PERMYY'), nc = numel(grdecl.PERMYY); end
   if nc == 0 && isfield(grdecl, 'PERMZZ'), nc = numel(grdecl.PERMZZ); end

   assert (nc > 0);

   vals = zeros([nc, 1]);
   comp = ones([3, 3]);

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   % Extract first row [k_{xx}, k_{xy}, k_{xz}] from input data.
   %
   if isfield(grdecl, 'PERMX') || isfield(grdecl, 'PERMXX') || isfield(grdecl, 'PERMH')
      if isfield(grdecl, 'PERMX')
         vals   = [vals, grdecl.PERMX ];
      elseif isfield(grdecl, 'PERMH')
         vals   = [vals, grdecl.PERMH];
      else
         vals   = [vals, grdecl.PERMXX];
      end
      comp(1,1) = size(vals,2);
      comp      = setDiagonalComponents(comp, 1, 2, 3);
   end

   if isfield(grdecl, 'PERMXY')
      vals      = [vals, grdecl.PERMXY];
      comp(1,2) = size(vals,2); comp(2,1) = comp(1,2);
   end

   if isfield(grdecl, 'PERMXZ')
      vals      = [vals, grdecl.PERMXZ];
      comp(1,3) = size(vals,2); comp(3,1) = comp(1,3);
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   % Extract second row [k_{yx}, k_{yy}, k_{yz}] from input data.
   %
   if isfield(grdecl, 'PERMYX')
      vals      = [vals, grdecl.PERMYX];
      comp(2,1) = size(vals,2); comp(1,2) = comp(2,1);
   end

   if isfield(grdecl, 'PERMY') || isfield(grdecl, 'PERMYY') || isfield(grdecl, 'PERMH')
      if isfield(grdecl, 'PERMY')
         vals   = [vals, grdecl.PERMY ];
      elseif isfield(grdecl, 'PERMH')
         vals   = [vals, grdecl.PERMH];
      else
         vals   = [vals, grdecl.PERMYY];
      end
      comp(2,2) = size(vals,2);
      comp      = setDiagonalComponents(comp, 2, 3, 1);
   end

   if isfield(grdecl, 'PERMYZ')
      vals      = [vals, grdecl.PERMYZ];
      comp(2,3) = size(vals,2); comp(3,2) = comp(2,3);
   end

   %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   % Extract third row [k_{zx}, k_{zy}, k_{zz}] from input data.
   %
   if isfield(grdecl, 'PERMZX')
      vals      = [vals, grdecl.PERMZX];
      comp(3,1) = size(vals,2); comp(1,3) = comp(3,1);
   end

   if isfield(grdecl, 'PERMZY')
      vals      = [vals, grdecl.PERMZY];
      comp(3,2) = size(vals,2); comp(2,3) = comp(3,2);
   end

   if isfield(grdecl, 'PERMZ') || isfield(grdecl, 'PERMZZ')
      if isfield(grdecl, 'PERMZ')
         vals   = [vals, grdecl.PERMZ ];
      else
         vals   = [vals, grdecl.PERMZZ];
      end
      comp(3,3) = size(vals,2);
      comp      = setDiagonalComponents(comp, 3, 1, 2);
   end
end

%--------------------------------------------------------------------------

function c = setDiagonalComponents(c, i, j, k)
   if c(j,j) == 1, c(j,j) = c(i,i); end
   if c(k,k) == 1, c(k,k) = c(i,i); end
end
