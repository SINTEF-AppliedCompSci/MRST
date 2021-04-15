function varargout = multiEig(varargin)
%Solve sequence of general (unsymmetric) eigenvalue problems using LAPACK
%
% SYNOPSIS:
%    d        = multiEig(A, sz)
%   [d, v]    = multiEig(A, sz)
%   [d, v, w] = multiEig(A, sz)
%
% PARAMETERS:
%   A  - Array (type `double`) containing the elements/entries of a sequence
%        of coefficient matrices--ordered consequtively.  Each matrix is
%        expected to be square and fairly small.
%
%   sz - Sequence of matrix block sizes.  The number of matrix blocks
%        contained in `A` is implicitly assumed to be `numel(rsz)`.
%
% RETURNS:
%   d - Array (type `double`, possibly complex) containing the eigenvalues of
%       the individual eigenvalue problems.
%
%   v - Array (type `double`) of right eigenvectors.  Optional.  Only returned
%       if requested.
%
%   w - Array (type `double`) of left eigenvectors.  Optional.  Only returned
%       if requested.
%
% NOTE:
%   The output in the single return value case is algorithmically equivalent
%   to the loop ::
%
%       pA = cumsum([1, sz .^ 2]);
%       pd = cumsum([1, sz]);
%       d  = zeros([pd(end) - 1, 1]);
%
%       for i = 1 : numel(sz),
%          iA = pA(i) : pA(i + 1) - 1;
%          id = pd(i) : pd(i + 1) - 1;
%
%          d(id) = eig(reshape(A(iA), [sz(i), sz(i)]), 'vector');
%       end
%
%   except for round-off errors.  Rank deficient matrices are supported.

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

if ispc,

   e = mexext('all');
   a = e(strcmp({ e.ext }, mexext)).arch;

   mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, 'microsoft', ...
               ['libmw', lib, '.lib']);

elseif isunix,

   mwlib = @(lib) ['-lmw', lib];

end

buildmex('-O', '-largeArrayDims', ...
         'multiEig.cpp',          ...
         mwlib('lapack'), mwlib('blas'));

% Now, run it.
[varargout{1:nargout}] = multiEig(varargin{:});
