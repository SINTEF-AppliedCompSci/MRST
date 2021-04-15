function varargout = interpRelPermTable(s, kr, si)
%Fast linear interpolation of tabulated relperm function (and derivative)
%
% SYNOPSIS:
%    kri        = interpRelPermTable(s, kr, si)
%   [kri, dkri] = interpRelPermTable(s, kr, si)
%
% PARAMETERS:
%   s    - Saturation nodes at which underlying function kr=kr(s) is
%          sampled.  The node values must be monotonically increasing.
%
%   kr   - Values of the underlying function kr=kr(s) at the nodes, s.
%          Must be level or increasing.
%
%   si   - Evaluation points (saturations) for new, interpolated, relperm
%          function values and derivatives.
%
% RETURNS:
%   kri  - Column vector of piecewise linearly interpolated values of the
%          relperm function kr=kr(s) at the evaluation points 'si'.
%
%   dkri - Column vector of piecewise constant derivative values of the
%          function kr=kr(s) at the evaluation points 'si'.
%          OPTIONAL.  Only computed and returned if specifically requested.
%
%          The derivative values are defined by the piecewise formulae
%
%                      '  0 ,      0   <= si <  s( 1 ),  % if 0 < s(1)
%               dkri = |  dm,   s( m ) <= si <  s(m+1),
%                      `  0 ,   s(end) <= si <=    1  ,
%
%           where dm = (kr(m+1) - kr(m)) / (s(m+1) - s(m)).
%
%           The caller needs to pay particular attention to the handling of
%           sub-interval end points when computing derivatives, as the
%           derivative is generally discontinuous.
%
% NOTE:
%   Function 'interpRelPermTable' is designed for efficient evaluation of
%   tabulated relative permeability curves defined by SWOF, SGOF, SOF2 and
%   related ECLIPSE keywords.
%
%   Using this utility function as a general interpolator is not advised.
%
% SEE ALSO:
%   `interpTable`, `dinterpTable`, `interp1q`, `dinterpq1`.

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


   s  = reshape(s , [], 1);   assert (all (diff(s) > 0));
   kr = reshape(kr, [], 1);
   si = reshape(si, [], 1);

   % Bin saturation values for efficient piecewise interpolation.
   %
   % Computational complexity:
   %
   %   \Omega(NUMEL(si) * LOG(NUMEL(s) + 2))) \approx \Omega(NUMEL(si))
   %
   [b,b] = histc(si, [-inf; s; inf]);    %#ok  % Ignore MLINT in M >= 7.9

   % Pad tables to cover saturation range [0,1].
   % Relperm end points equal to min/max kr values, respectively.
   %
   S  =    [s(1)-1;   s  ; s(end)+1] ;
   Kr = kr([   1  , 1:end,   end   ]);

   % Compute piecewise linear interpolant coefficients...
   DS = diff(S);
   T  = (si - S(b)) ./ DS(b);

   % ... and derive interpolated relperm values.
   varargout{1} = T.*Kr(b+1) + (1-T).*Kr(b);

   if nargout > 1,
      % Caller requested derivatives too.
      DKr = diff(Kr) ./ DS;
      varargout{2} = DKr(b);
   end
end
