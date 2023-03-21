function zcorn = checkAndRepairZCORN(zcorn, cartDims, varargin)
%Detect and repair artifacts that may occur in corner-point specification.
%
% SYNOPSIS:
%   zcorn = checkAndRepairZCORN(zcorn, dims)
%   zcorn = checkAndRepairZCORN(zcorn, dims, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function `checkAndRepairZCORN` detects and repairs the following, rare,
%   conditions
%
%      - Upper corners of a cell below lower corners of that same cell
%      - Lower corners of a cell below that cell's lower neighbour's upper
%        corners.
%
%   These will typically arise as a result of finite precision output from
%   a corner-point grid generator.
%
%   The repair strategy, if applicable, is as follows,
%
%     - If an upper corner is below the corresponding lower corner on the
%       same pillar, then the upper corner depth is assigned to be equal to
%       the lower corner depth.
%
%     - If a cell's lower corner is below that cell's lower neighbour's
%       upper corner on the same pillar, then the lower neighbour's upper
%       corner is assigned the corner depth of the upper cell's lower
%       corner.
%
% PARAMETERS:
%   zcorn   - Corner-point depth specification of a corner-point grid.
%             This value typically corresponds to the 'ZCORN' field of a
%             data structure created by function 'readGRDECL'.
%
%   dims    - Cartesian dimensions of the corner-point geometry.  Assumed
%             to be a three-element vector of (positive) extents.
%             Typically corresponds to the field 'cartDims' of a 'grdecl'
%             structure created by function 'readGRDECL'.
%
% KEYWORD ARGUMENTS:
%
%   Active  - A `prod(dims)`-element vector of active cell mappings.
%              Zero/false signifies inactive cells while non-zero/true
%              signifies active cells. Typically corresponds to the field
%              `ACTNUM` of a 'grdecl' structure.
%                          Default value: `Active=[]` (-> All cells active).
%
%   Verbose - Whether or not to emit informational messages during the
%             computational process. 
%             Default value: `Verbose = mrstVerbose()`.
%
% RETURNS:
%   zcorn   - Corner-point depth specification for which identified
%             problems have been corrected.  If there are no problems, then
%             this is the same as input array 'zcorn'.
%
% NOTE:
%   This function is used to implement option `RepairZCORN` of the
%   corner-point processor, `processGRDECL`.
%
% SEE ALSO:
%   `readGRDECL`, `processGRDECL`, `mrstVerbose`

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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


   assert (isnumeric(zcorn) && isnumeric(cartDims), ...
           'Input parameters ''zcorn'' and ''cartDims'' must be numeric.');

   assert (numel(cartDims) == 3, ...
          ['Function ''%s'' is only supported in three ', ...
           'space dimensions'], mfilename);

   assert (numel(zcorn) == prod(2 * cartDims), ...
          ['Corner point depths must be specified independently for\n', ...
           'each of eight corners in each of %dx%dx%d=%d global cells.', ...
           '\nExpected %d corner point depths but got %d.'], ...
           cartDims(1), cartDims(2), cartDims(3), prod(cartDims), ...
           prod(2 * cartDims), numel(zcorn));

   opt = struct('active', [], 'verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   act = get_active(opt, cartDims);
   Z   = reshape(zcorn, 2 * cartDims);

   [Z, modtbb]  = top_below_bottom(Z, act);
   [Z, modbblt] = bottom_below_lower_top(Z, act);

   if modtbb || modbblt,
      % Either (or both) of the conditions modified the CP depths.  Write
      % modified result back to original ZCORN array.
      zcorn(:) = Z;

      dispif(opt.verbose && modtbb, ...
             'Modified %d depths violating top-below-bottom\n', modtbb);

      dispif(opt.verbose && modbblt, ...
             'Modified %d depths violating bottom-below-lower-top\n', ...
             modbblt);
   end
end

%--------------------------------------------------------------------------

function act = get_active(opt, cartDims)
   act = opt.active;

   if isempty(act),

      act = true(cartDims);

   else
      if ~ (isnumeric(act) || islogical(act)),
         error('Active maps of type ''%s'' are not supported', class(act));
      end

      assert (numel(act) == prod(cartDims), ...
              'Active map must specify one indicator for each global cell');

      if isnumeric(act),

         act = reshape(act ~= 0, cartDims);

      else
         % ISLOGICAL(act)

         act = reshape(act, cartDims);

      end
   end

   % Produce [1, 1, 2, 2, 3, 3, ...,        ...
   %          cartDims(d)-1, cartDims(d)-1, ...
   %          cartDims(d)  , cartDims(d)]
   %
   % Equivalently:
   %   ind = @(d) reshape(repmat(1 : cartDims(d), [2, 1]), 1, [])
   %   ind = @(d) rldecode(1 : cartDims(d), 2, 2)
   %
   % Using 'fix' has the least computational overhead in this case.
   %
   ind = @(d) fix((2 : 2*cartDims(d) + 1) ./ 2);

   % Expand active map to cover all pillars repeated for all cells.
   act = act(ind(1), ind(2), :);
end

%--------------------------------------------------------------------------

% Check if any cell's upper corner-points are below that cell's lower
% corner-points on the same pillar.  Assign upper=lower in that case.
function [Z, modified]  = top_below_bottom(Z, act)
   tbb = Z(:, :, 1:2:end) > Z(:, :, 2:2:end);  % Top below bottom

   assert (all(size(tbb) == size(act)), 'tbb: Internal error');

   invalid = find(tbb & act);
   if ~ isempty(invalid),
      % There are some, although hopefully rather few, points of trouble.

      [I, J, K] = ind2sub(size(tbb), invalid);

      dst = sub2ind(size(Z), I, J, 2*K - 1);   % Upper corner point
      src = sub2ind(size(Z), I, J, 2*K    );   % Lower corner point

      Z(dst)   = Z(src);         % Assign upper=lower at troublesome points
      modified = numel(invalid); % Number of modified depths
   else
      modified = false;
   end
end

%--------------------------------------------------------------------------

% Check if any cell's lower/bottom corner-points are below that cell's
% lower neighbour's upper corner-points on the same pillar.  Assign
% "lower's upper = upper's lower" in that case.
function [Z, modified] = bottom_below_lower_top(Z, act)
   bblt = Z(:, :, 2:2:end-1) > Z(:, :, 3:2:end);  % Bottom below lower top

   act_face = act(:, :, 1:end-1) & act(:, :, 2:end);
   invalid  = find(bblt & act_face);  % BBLT & upper/lower both active
   if ~ isempty(invalid),
      % There are some, although hopefully rather few, points of trouble.

      [I, J, K] = ind2sub(size(bblt), invalid);

      dst = sub2ind(size(Z), I, J, 2*K + 1);   % Lower's upper
      src = sub2ind(size(Z), I, J, 2*K    );   % Upper's lower

      Z(dst)   = Z(src);         % Assign lower's upper=upper's lower
      modified = numel(invalid); % Number of modified depths
   else
      modified = false;
   end
end
