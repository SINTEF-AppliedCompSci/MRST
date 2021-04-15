function [yi, dyidx] = interpTable(X, Y, xi, varargin)
%Interpolate a one-dimensional table, possibly using splines.
%
% SYNOPSIS:
%   yi = interpTable(X, Y, xi)
%   yi = interpTable(X, Y, xi, 'pn1', pv1, ...)
%
% PARAMETERS:
%   X       - Nodes at which underlying function y=y(x) is sampled.
%
%   Y       - Values of the underlying function y=y(x) at the nodes, `X`.
%
%   xi      - Evaluation points for new, interpolated, function values.
%
% OPTIONAL PARAMETERS:
%   'spline' - Whether or not to use spline interpolation.
%              Logical.  Default value: `spline=false` (use linear
%              interpolation/extrapolation).
%
% RETURNS:
%   yi - Approximate (interpolated/extrapolated) values of the function
%        y=y(x) at the points `xi`.
%
% SEE ALSO:
%   `dinterpTable`, `interp1`, `spline`.

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


   opt = struct('spline', false);
   opt = merge_options(opt, varargin{:});

   if opt.spline
      method = 'spline';
   else
      method = 'linear';
   end

   yi = interp1(X, Y, xi, method, 'extrap');
   if nargout > 1
       dyidx = dinterpTable(X, Y, xi, varargin{:});
   end
end
