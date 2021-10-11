function [kr, pc] = tabulatedSatFunc(T)
%Construct rel-perm and capillary pressure evaluators from tabulated data
%
% SYNOPSIS:
%   [kr, pc] = tabulatedSatFunc(table)
%
% PARAMETERS:
%   table - Tabulated data for relative permeability and capillary pressure
%           as a function of water saturation.  The table must be a three-
%           or four-column array in which the first column is interpreted
%           as water saturation, the second is water relative permeability
%           and the third is relative permeability of oil.  The fourth
%           column, if present, is interpreted as the oil-water capillary
%           pressure function, P_oil-P_wat.
%
%           The first column, saturation, must be monotonically increasing.
%           The second column must be level or increasing, while the third
%           column must be level or decreasing.  The fourth column, if
%           present, must be level or decreasing.
%
% RETURNS:
%   kr - Function handle that supports the following calling convention
%
%            KR       = kr(s)
%           [KR, dKR] = kr(s)
%
%        in which the first syntax evaluates the relative permeability
%        function, using linear interpolation, in 'table' at the saturation
%        point 's' and the second syntax additionally provides derivatives
%        of the rel-perm function (piecewise constant interpolation).
%
%        The return value 'KR' has one row (and two columns) for each *row*
%        in the saturation 's'.  The return value 'dKR', if requested, has
%        one row and four colums for each row in the saturation 's'.  The
%        column values are the partial derivatives
%
%             [  d KR(:,1) / d s(:,1)  ,  ...
%                d KR(:,2) / d s(:,1)  ,  ...
%                d KR(:,1) / d s(:,2)  ,  ...
%                d KR(:,2) / d s(:,2)  ]
%
%        respectively.  The diagonal of the Jacobian matrix at each
%        saturation value may thus be extracted as 'dKR(:, [1, 4])'.
%
%        For compatiblity with the rest of MRST, the 'kr' function accepts
%        any positive number of input arguments.  All arguments other than
%        the first are ignored.
%
%   pc - Function handle that supports the following calling convention
%
%            PC       = pc(s)
%           [PC, dPC] = pc(s)
%
%        The return value 'PC' is the capillary pressure P_oil - P_wat
%        evaluated at the saturation point 's'.  The second return value,
%        if requested, is the derivative of capillary pressure with respect
%        to s(:,1).
%
% NOTE:
%   The evaluators constructed by function 'tabulatedSatFunc' are
%   explicitly two-phase only.  Specifically, the evaluators use only the
%   first column of the input parameter 's' and explicitly assume that the
%   second column, if present, satisfies the saturation constraint
%
%             s(:,2) = 1 - s(:,1)

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


   check_table(T);

   kr = @(s, varargin) interp_relperm(s, T(:, 1 : 3));

   if size(T, 2) == 4,
      pc = @(s, varargin) interp_pcap(s, T(:, [1, 4]));
   else
      pc = @zero_capillary;
   end
end

%--------------------------------------------------------------------------

function check_table(T)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:3) < 0)), ...
           'Saturation and rel-perm values must be non-negative.');

   assert (~ (T( 1 ,2) > 0), ...
           'Water must be immobile at connate water saturation.');
   assert (~ (T(end,3) > 0), ...
           'Oil must be immobile at maximum water saturation.');

   assert (all (diff(T(:,1)) > 0), ...
           'Water saturation must be monotonically increasing.');
   assert (~any(diff(T(:,2)) < 0), ...
           'Water rel-perm must be level or increasing down column.');
   assert (~any(diff(T(:,3)) > 0), ...
           'Oil rel-perm must be level or decreasing down column.');
end

%--------------------------------------------------------------------------

function varargout = interp_relperm(s, T)
   sw = max(s(:,1), T(1,1));

   [varargout{1 : nargout}] = interpolate(sw, T(:, 1), T(:, 2 : 3));

   if nargout > 1,
      varargout{2} = [ varargout{2}(:, 1)    , ...
                                               ...
                      zeros([size(sw, 1), 2]), ...
                                               ...
                      -varargout{2}(:, 2) ];
   end
end

%--------------------------------------------------------------------------

function varargout = interp_pcap(s, T)
   [varargout{1 : nargout}] = interpolate(max(s(:,1), T(1,1)), ...
                                          T(:, 1), T(:, 2));
end

%--------------------------------------------------------------------------

function varargout = zero_capillary(s, varargin)
   varargout(1 : nargout) = { zeros([size(s, 1), 1]) };
end

%--------------------------------------------------------------------------

function varargout = interpolate(sw, s, y)
   % Bin saturation values for efficient piecewise interpolation.
   %
   % Computational complexity:
   %
   %   \Omega(NUMEL(sw) * LOG(NUMEL(s) + 2))) \approx \Omega(NUMEL(sw))
   %
   [b, b] = histc(sw, [-inf; s(2 : end-1); inf]);                      %#ok

   % Compute piecewise linear interpolant coefficients...
   DS = diff(s);
   T  = (sw - s(b)) ./ DS(b);

   % ... and derive interpolated values.
   varargout{1} = bsxfun(@times,   T  , y(b + 1, :)) + ...
                  bsxfun(@times, 1 - T, y(  b  , :));

   if nargout > 1,
      % Caller requested derivatives too.
      DY = bsxfun(@rdivide, diff(y), DS);
      varargout{2} = DY(b, :);
   end
end
