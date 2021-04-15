function varargout = blockDiagIndex(m, n)
%Compute subscript or linear index to nonzeros of block-diagonal matrix
%
% SYNOPSIS:
%   [i, j] = blockDiagIndex(m, n)
%   [i, j] = blockDiagIndex(m)
%    i     = blockDiagIndex(...)
%
% PARAMETERS:
%   m,n    - Vectors of dimensions/block sizes of each diagonal block.
%            Input `m` is interpreted as the row dimension of each block
%            while `n` is interpreted as the column dimension of each
%            block.
%
%            Vector `n` must have the same number of entries as vector `m`,
%            i.e., the number of blocks.
%
%            If `n` is not supplied, function `blockDiagIndex` will behave
%            as if it were called as::
%
%                [i, j] = blockDiagIndex(m, m)
%
%            In other words, using square blocks.
%
% RETURNS:
%   i,j    - Index vectors that may be used to form sparse block-diagonal
%            matrix.  This process is demonstrated in the example below.
%
%            The return values are of type `double` irrespective of the
%            class of the dimension vectors `m` and `n`.
%
%            If called using a single output value, that value is the
%            linear index, computed using `sub2ind`, of the index pair
%            (i,j).
%
% EXAMPLE:
%   m      = [ 1 ; 2 ; 3 ];
%   n      = [ 2 ; 3 ; 4 ];
%   [i, j] = blockDiagIndex(m, n);
%   A      = sparse(i, j, 1 : sum(m .* n));
%   full(A), spy(A)
%
% SEE ALSO:
%   `rldecode`, `mcolon`, `sub2ind`.

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

   mrstNargInCheck(1, 2, nargin);

   if nargin == 1, n = m; end

   assert (numel(n) == numel(m), ...
           'Dimension vectors must have the same number of entries.');

   m = reshape(double(m), [], 1);
   n = reshape(double(n), [], 1);

   pos = cumsum([1 ; m]);
   p1  = pos(1 : end - 1);
   p2  = pos(2 : end) - 1;
   i   = mcolon(rldecode(p1, n),  rldecode(p2, n)) .';
   j   = rldecode((1 : sum(n)) .', rldecode(m, n));

   if nargout == 1,
      varargout{1} = sub2ind([sum(m), sum(n)], i, j);
   else
      varargout{1} = i;
      varargout{2} = j;
   end
end
