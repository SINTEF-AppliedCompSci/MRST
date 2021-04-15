function varargout = invv(varargin)
%Compute inverse of sequence of square matrices using LAPACK.
%
% SYNOPSIS:
%   B = invv(A, sz)
%
% PARAMETERS:
%   A  - Column array (type `double`) containing the elements/entries of a
%        sequence of square matrices--ordered consequtively.  Each matrix
%        is expected to be small (dimension less than 50).
%
%   sz - Sequence of matrix block sizes, i.e. the number of rows/columns in
%        each square matrix. The number of matrix blocks contained in `A`
%        is implicitly assumed to be `numel(sz)`.
%
% RETURNS:
%   B  - Array (type `double`) containing the elements/entries of the
%        sequence of inverses of the square matrices contained in `A` and
%        whose block sizes are `sz`.
%
% NOTE:
%   The output is algorithmically equivalent to the loop ::
%
%       p = cumsum([1, sz.^2]);
%       B = zeros([p(end) - 1, 1]);
%       for i = 1 : numel(sz),
%          ix    = p(i) : p(i + 1) - 1;
%          B(ix) = inv(reshape(A(ix), sz([i, i])));
%       end
%
%   except for round-off errors.
%
% EXAMPLE:
%   a = [1; 0; ... % First row, first matrix
%        0; 1; ... % Second row, first matrix
%        2; 3; ... % First row, second matrix
%        4; 5];    % Second row, second matrix
%   b = [2, 2];
%   result = invv(a, b)

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


if ispc

   e = mexext('all');
   a = e(strcmp({ e.ext }, mexext)).arch;

   mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, 'microsoft', ...
               ['libmw', lib, '.lib']);

elseif isunix

   mwlib = @(lib) ['-lmw', lib];

end

buildmex('invv.c', mwlib('lapack'), mwlib('blas'));

% Now, run it.
[varargout{1:nargout}] = invv(varargin{:});
