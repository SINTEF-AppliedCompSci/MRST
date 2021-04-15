function [A, n] = rlencode(A, dim)
%Compute run length encoding of array A along dimension dim.
%
% SYNOPSIS:
%   [A, n] = rlencode(A)
%   [A, n] = rlencode(A, dim)
%
% PARAMETERS:
%   A         - Array.  Numeric or cell array of string ("cellstring").
%   dim       - dimension of `A` where run length encoding is done.
%               `dim > 0`, `dim <= ndims(A)`.
%               OPTIONAL.  Default value: `dim=1`.
%
% RETURNS:
%   A         - Compressed `A` where repeated layers are removed.
%   n         - repetition count of repeated layers in original `A`.
%
% EXAMPLE:
%   % 1) Regular numeric matrix
%   A      = [ 1, 2, 3, 4 ; ...
%              1, 2, 3, 4 ; ...
%              3, 4, 5, 6 ; ...
%              3, 3, 3, 3 ; ...
%              3, 3, 4, 5 ; ...
%              3, 3, 4, 5 ];
%
%   [A, n] = rlencode(A, 1);
%
%   assert (isequal(A, [ 1, 2, 3, 4 ; ...
%                        3, 4, 5, 6 ; ...
%                        3, 3, 3, 3 ; ...
%                        3, 3, 4, 5 ]))
%
%   assert (isequal(n, [2, 1, 1, 2] .'))
%
%   % 2) Cell array of string
%   S = { 'aaaa' ; 'aaaa' ; 'bb' ; 'cccc' ; 'cccc' ; 'd'; 'd'; 'd' ; 'd' };
%   [s, n] = rlencode(S);
%
%   assert (isequal(s, { 'aaaa' ; 'bb'; 'cccc'; 'd' }))
%   assert (isequal(n, [ 2; 1; 2; 4 ]))
%
% SEE ALSO:
%   `rldecode`.

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

   if isempty(A)
      n = 0;
      return
   end

   if nargin < 2
      dim = 1;
   end

   % Pick out (1:end - 1, :, ...) and (2:end, :, ...) in a multidimensional
   % way.  In particular, swap dimensions 1 and 'dim' so that we only have
   % to consider compressing along the first dimension.
   d = 1 : ndims(A) ;   d([1, dim]) = [dim, 1];
   B = permute(A, d);

   % Find positions where layers differ.
   i = [find(any(different_layers(B), 2)) ; size(B, 1)];

   % Compare differences in position to find run length.
   n = diff([0; i]);

   % Re-swap dimensions 1 and 'dim' to form return value.
   sz = size(B); sz(1) = numel(i);
   A  = permute(reshape(B(i,:), sz), d);
end

%--------------------------------------------------------------------------

function d = different_layers(B)
   top = B(1 : (end - 1), :);
   bot = B(2 : (end - 0), :);

   if isnumeric(B)
      d = xor(top ~= bot, isnan(top) & isnan(bot));
   else
      % Handles case of cellstring and, in MATLAB >= R2016b, STRING too.
      d = ~strcmp(top, bot);
   end
end
