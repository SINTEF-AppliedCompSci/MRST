function varargout = matrixBlocksFromSparse(varargin)
%Extract block-diagonal matrix elements from sparse matrix
%
% SYNOPSIS:
%   Sbd = matrixBlocksFromSparse(S, sz)
%
% PARAMETERS
%   S  - Sparse matrix.  Assumed to be square and block-diagonal with
%        square blocks.  The diagonal blocks may contain zero elements.
%
%   sz - Block sizes.  An n-by-1 INT32 vector of positive numbers such that
%        the i'th diagonal block of 'S' is sz(i)-by-sz(i).
%
% RETURNS:
%   Sbd - Block-diagonal matrix elements.  A SUM(sz.^2)-by-1 array of real
%         numbers.  May contain zeros.
%
% NOTE:
%   If each individual block along the diagonal is full, then the output is
%   equivalent to nonzeros(S).
%
% SEE ALSO:
%   `nonzeros`, `computeMultiPointTrans`, `incompMPFA`.

%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

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


   % If we end up here, matrixBlocksFromSparse needs to be built.
   % Next time, this will not happen.
   buildmex -largeArrayDims matrixBlocksFromSparse.c

   % Now, run it.
   [varargout{1:nargout}] = matrixBlocksFromSparse(varargin{:});
end
