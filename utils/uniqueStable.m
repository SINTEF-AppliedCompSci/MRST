function [c, ia, ic] = uniqueStable(a, varargin)
%Support `unique(A, 'stable')` in all versions of MATLAB
%
% SYNOPSIS:
%   [c, ia, ic] = uniqueStable(a)
%   [c, ia, ic] = uniqueStable(a, 'rows')
%   [c, ia, ic] = uniqueStable(..., 'use_fallback')
%
% DESCRIPTION:
%   This is a pure MATLAB compatibility implementation of ::
%
%        unique(A,         'stable')
%        unique(A, 'rows', 'stable')
%
%   for MATLABs prior to release R2012a (MATLAB 7.14).  In versions 7.14
%   and later, this simply forwards the parameters to the built-in version
%   of function `unique`.
%
%   Function `unique`'s 'stable' option returns the unique values in the
%   input in the same order as they appear in the input.  Without 'stable',
%   the unique elements are returned in sorted order.
%
% PARAMETERS:
%   a      - Numeric array.  The fall-back implementation does not support
%            cellstrings.
%
%   'rows' - Exact string 'rows' indicating that we should compute unique
%            rows of 'a' rather than merely unique single elements.
%
%   'use_fallback' -
%            Exact string 'use_fallback'.  This is mainly intended for
%            testing and development purposes.  The option bypasses the
%            MATLAB version check logic and forces the use of the fall-back
%            implementation.
%
% RETURNS:
%   c  - Unique elements (or rows) from input array 'a', in order of
%        appearance in the input array.
%
%   ia - Index into input 'a' such that ALL(ALL(c == a(ia, :)))
%
%   ic - Index into output 'c' such that ALL(ALL(c(ic, :) == a)).
%
% NOTE:
%   This function uses `sortrows`.
%
% EXAMPLE:
%   % 1) Unique elements in order of appearance
%   a = [ 3, 3, 3 ; 1, 4, 1 ; 2, 2, 2 ; 1, 1, 1 ; 2, 2, 2 ; 3, 3, 3 ];
%
%   [c, ia, ic] = uniqueStable(a);
%   assert (all(c == [ 3 ; 1 ; 2 ; 4 ]), ...
%           'Stable sort regression in ''uniqueStable''.')
%   assert (all(c == a(ia)), '''IA'' Regression in ''uniqueStable''.')
%   assert (all(all(reshape(c(ic), size(a)) == a)), ...
%           '''IC'' Regression in ''uniqueStable''.')
%
%   % 2) Unique matrix rows in order of appearance
%   [c, ia, ic] = uniqueStable(a, 'rows');
%   assert (all(c(:,1) == [ 3 ; 1 ; 2 ; 1 ]), ...
%           'Stable sort regression in ''uniqueStable/rows''.')
%   assert (c(2,2) == 4, ...
%           'Stable sort regression ((2,2)==4) in ''uniqueStable/rows''.')
%   assert (all(all(c == a(ia, :))), ...
%           '''IA'' Regression in ''uniqueStable/rows''.')
%   assert (all(all(c(ic, :) == a)), ...
%           '''IC'' Regression in ''uniqueStable/rows''.')
%
% SEE ALSO:
%   `unique`, `sortrows`.

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
   persistent old_version
   if isempty(old_version)
       if mrstPlatform('matlab')
           old_version = ~exist('verLessThan', 'file') || verLessThan('matlab', '7.14');
       else
           % Assume recent octave for the time being.
           old_version = false;
       end
   end
   if isempty(a)
      [c, ia, ic] = deal([]);
      return
   end

   if force_fallback(varargin{:}) || old_version

      % Caller explicitly requested that the fall-back option be used or
      % we're currently targeting an older release of MATLAB.  Feature
      % UNIQUE(..., 'stable') was introduced in R2012a (a.k.a MATLAB 7.14).

      [c, ia, ic] = fall_back(a, varargin{:});

   else
      % MATLAB is sufficiently recent.  Just invoke built-in function.

      [c, ia, ic] = unique(a, varargin{:}, 'stable');

   end
end

%--------------------------------------------------------------------------

function b = force_fallback(varargin)
   b = (nargin > 0) && iscellstr(varargin) && ...
      any(strcmpi(varargin, 'use_fallback'));
end

%--------------------------------------------------------------------------

function [c, ia, ic] = fall_back(a, varargin)
   opt   = get_options(varargin{:});
   isrow = is_row_shape(a);

   if isrow && opt.rows
      % UNIQUE(A, 'rows', 'stable') on row vector implies
      %
      %    C = A, IA = 1, IC = 1

      [c, ia, ic] = deal(a, 1, 1);

   else

      if opt.rows && (ndims(a) > 2)                             %#ok<ISMAT>
         error('uniqueStable:Rows:NDArray', ...
               'Option ''rows'' is only supported for matrices.');
      end

      if isrow || ~opt.rows
         % Row vector or 2D array without 'rows' option.  Linearise to
         % treat as column vector.
         a = reshape(a, [], 1);
      end

      t      = sortrows([a, (1 : size(a, 1)) .']);  % Stable sort.
      [n, n] = rlencode(t(:, 1 : (end - 1)), 1);                %#ok<ASGLU>
      pos    = cumsum([ 1 ; n ]);

      % Extract unique values and reorder according to original position.
      x      = t(pos(1 : end - 1), :);
      [i, i] = sort(x(:, end));                                 %#ok<ASGLU>

      c  = x(i, 1 : (end - 1));
      ia = x(i, end);

      renum         = zeros(size(ia));
      renum(i)      = 1 : numel(renum);
      ic            = zeros([size(t, 1), 1]);
      ic(t(:, end)) = rldecode(renum, n);

      if isrow
         % Input was row vector.  Return row shape.
         c = reshape(c, 1, []);
      end
   end
end

%--------------------------------------------------------------------------

function opt = get_options(varargin)
   opt = struct('rows', false);

   if (nargin > 0) && iscellstr(varargin) && ...
         any(strcmp(varargin, 'rows'))
      opt.rows = true;
   end
end

%--------------------------------------------------------------------------

function row_p = is_row_shape(a)
   sz = size(a);

   row_p = (numel(sz) == 2) && (sz(1) == 1);
end
