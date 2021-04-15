function varargout = matrixBlocksFromSparse(varargin)
%Extract block-diagonal matrix elements from sparse matrix
%
% SYNOPSIS:
%    Sbd       = matrixBlocksFromSparse(S, rsz)
%    Sbd       = matrixBlocksFromSparse(S, rsz, csz)
%   [Sbd, pos] = matrixBlocksFromSparse(...)
%
% PARAMETERS
%   S   - Sparse matrix.  Assumed to be block-diagonal.  The diagonal blocks
%         may contain zero elements.
%
%   rsz - Row sizes.  An n-by-1 vector of positive integers such that the
%         number of rows of block 'i' is 'rsz(i)'.
%
%   csz - Column sizes.  An n-by-1 vector of positive integers such that the
%         number of columns of block 'i' is 'csz(i)'.
%
%         Optional.  If unspecified, copied from 'rsz'--in other words, if
%         'csz' is not specified, the return value will be constructed from
%         square blocks.
%
% RETURNS:
%   Sbd - Block-diagonal matrix elements.  A SUM(rsz .* csz)-by-1 array of
%         real numbers.  May contain zeros.
%
%   pos - Start pointers for individual blocks.  In particular, the 'i'-th
%         block of 'S' is contained in entries
%
%            pos(i) : pos(i + 1) - 1
%
%         of 'Sbd'.  This is a convenience value only as the same pointer
%         structure can be computed using the statement
%
%            pos = cumsum([1 ; rsz .* csz])
%
% NOTE:
%   If each individual block along the diagonal is full, then the output is
%   equivalent to nonzeros(S).
%
% SEE ALSO:
%   `nonzeros`, `computeMultiPointTrans`, `incompMPFA`.

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


   % If we end up here, matrixBlocksFromSparse needs to be built.
   % Next time, this will not happen.
   buildmex -O -largeArrayDims matrixBlocksFromSparse.c

   % Now, run it.
   [varargout{1:nargout}] = matrixBlocksFromSparse(varargin{:});
end
