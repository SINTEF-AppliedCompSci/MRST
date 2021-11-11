function x = mcolon(lo, hi, s)
%Compute vector of consecutive indices from separate lower/upper bounds.
%
% SYNOPSIS:
%   ind = mcolon(lo, hi)
%   ind = mcolon(lo, hi, s)
%
% PARAMETERS:
%   lo  - Vector of start values (lower bounds).
%   hi  - Vector of end values (upper bounds).
%   s   - Vector of strides.
%         Optional.  Default value: s = ones([numel(lo), 1]) (unit stride).
%
% RETURNS:
%   ind - `[lo(1):hi(1)     , lo(2):hi(2)     ,..., lo(end):hi(end)]`
%   ind - `[lo(1):s(1):hi(1), lo(2):s(2):hi(2),...,lo(end):s(end):hi(end)]`
%
% EXAMPLE:
%   lo  = [1 1 1 1]; hi = [2 3 4 5];
%   ind = mcolon(lo, hi)
%
% NOTES:
%   Note that `ind` has type `double` irrespective of the type of its input
%   parameters.
%
%   `mcolon` may be implemented in terms of `arrayfun` and `horzcat`, e.g., ::
%
%      ind = arrayfun(@colon, lo, hi, 'UniformOutput', false);
%      ind = [ind{:}];
%   or ::
%
%      ind = arrayfun(@colon, lo, s, hi, 'UniformOutput', false);
%      ind = [ind{:}];
%   but the current implementation is faster.

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

   % Support general input class.
   lo = double(lo);
   hi = double(hi);

   if numel(lo) == 1, lo = repmat(lo, size(hi)); end
   if numel(hi) == 1, hi = repmat(hi, size(lo)); end

   if numel(lo) ~= numel(hi)
      error('Dimension:Mismatch', ...
            'LO and HI must have the same number of elements');

   elseif isempty(lo)
      % Empty index boundaries -> empty result.
      x = [];

   elseif nargin < 3
      % x = mcolon(lo, hi), stride = 1.
      x = unit_stride(lo, hi);

   else
      % x = mcolon(lo, hi, s), stride = s.
      x = user_defined_stride(lo, hi, s);
   end
end

%--------------------------------------------------------------------------

function x = unit_stride(lo, hi)
   % Remove lo-hi pairs where numel(lo:hi)==0
   i    = hi >= lo;

   if ~ any(i), x = []; return, end

   hi   = hi(i);
   lo   = lo(i);
   m    = numel(lo);           % Number of bins
   d    = double(hi - lo + 1); % Number of elements per bin (hi/lo included)
   n    = sum(d);              % Total number of (expanded) elements

   % Preallocate result.  Recall: CUMSUM(ONES([1,n])) == 1:n
   x    = ones([1, n]);
   x(1) = lo(1); % Starting index of first bin in first position.

   % Prepare for running sum accumulation.  The following properties are
   % used
   %
   %   1) 1 + cumsum(d(1:end-1)) is starting positions of bins 2:m
   %   2) lo(2:m) - hi(1:m-1) is offset from hi(1:m-1) to reach lo(2:m)
   %
   x(1 + cumsum(d(1:end-1))) = lo(2:m) - hi(1:m-1);

   % Alternative setups to previous two statements.
   %   x(cumsum([1 , d(1:end-1)])) = [ lo(1) , lo(2:m) - hi(1:m-1) ]
   %   x(cumsum([1 , d(1:end-1)])) = lo - [ 0 , hi(1:m-1) ]

   % Compute result by running sum of index accumulators (mostly ONES).
   %
   % Note: Running sum is adjusted (i.e., reset to lo(i)) at each bin
   % boundary 'i' by 2) above.
   x = cumsum(x);
end

%--------------------------------------------------------------------------

function x = user_defined_stride(lo, hi, s)
   s = double(s);
   if numel(s) == 1
      s = repmat(s, size(lo));
   end

   % Remove lo-hi-s triplets where numel(lo:s:hi)==0
   i    = ((hi >= lo) & (s > 0)) | ((hi <= lo) & (s < 0));

   if sum(i) == 0, x = []; return, end

   [hi, lo, s] = deal(hi(i), lo(i), s(i));

   % Compute lo + (0:(hi-lo)/stride)*stride
   % Fix or hack: avoid roundoff error in floor((hi-lo)./s) when hi-lo = N*s
   % for some natural number N.
   e    =  (1 - 2*(hi < lo)) * eps;
   hi   = fix((e + hi - lo) ./ s);

   m    = numel(lo);
   d    = double(hi + 1);
   n    = sum(d);

   assert (all(d > 0), ...
          ['Internal error in ''%s'': Bins with non-positive ', ...
           'number of elements detected'], mfilename);

   ind = 1 + cumsum(d(1:end-1));

   % Expand lo to [lo(1) lo(1) ... lo(2) lo(2) ... lo(end)]
   LO      = zeros(1, n);
   LO(1)   = lo(1);
   LO(ind) = lo(2:m) - lo(1:m-1);
   LO      = cumsum(LO);

   % Expand stride
   S      = zeros(1, n);
   S(1)   = s(1);
   S(ind) = s(2:m) - s(1:m-1);
   S      = cumsum(S);

   x      = ones(1, n);
   x(1)   = 0;
   x(ind) = -hi(1:m-1);

   x = cumsum(x).*S + LO;
end
