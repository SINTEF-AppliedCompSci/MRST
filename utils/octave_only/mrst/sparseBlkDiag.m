function A = sparseBlkDiag(varargin)
%Form sparse block diagonal matrix from individual blocks.
%
% SYNOPSIS:
%   A = sparseBlkDiag(A11, A22, ..., Ann)
%
% PARAMETERS:
%   A11, A22, ..., Ann --
%       Matrix blocks, sparse or otherwise, which will be diagonally
%       concatenated to produce the output matrix.
%
% RETURNS:
%   A = | A11  0   ...  0  |
%       |  0  A22       0  |
%       |  :   :    \   :  |
%       |  0   0   ... Ann |
%
%   The return value, A, is a sparse matrix irrespective of the sparseness
%   of the input blocks A11, A22, ..., Ann.
%
% EXAMPLE:
%   A = [2, 1; 1, 2];
%   A = sparseBlkDiag(A, A, A, A);
%   full(A)
%
% NOTE:
%   This implementation is based on recursion.  The function might not be
%   appropriate in problems where there are significant memory constraints,
%   even considering copy-on-write optimizations behind the scenes.
%
% SEE ALSO:
%   blkdiag.

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


   if     nargin == 0,
      A = sparse(0, 0);
      return;
   elseif nargin == 1,
      A11 = varargin{1};
      A22 = sparse(0, 0);
   elseif nargin == 2,
      A11 = varargin{1};
      A22 = varargin{2};
   else
      p   = fix((1 + nargin) / 2);
      A11 = sparseBlkDiag(varargin{0 + 1 :  p });
      A22 = sparseBlkDiag(varargin{p + 1 : end});
   end

   if isempty(A11), A11 = sparse(0,0); end
   if isempty(A22), A22 = sparse(0,0); end

   if ~issparse(A11), A11 = sparse(A11); end
   if ~issparse(A22), A22 = sparse(A22); end

   [m1, n1] = size(A11);
   [m2, n2] = size(A22);

   A = [     A11      , sparse(m1, n2); ...
        sparse(m2, n1),      A22     ];
end
