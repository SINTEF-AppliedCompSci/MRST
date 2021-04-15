function varargout = leastsq_svd(varargin)
%Solve sequence of linear least squares problems using SVD method
%
% SYNOPSIS:
%   x = leastsq_svd(A, b, rsz, csz)
%
% PARAMETERS:
%   A   - Array (type `double`) containing the elements/entries of a sequence
%         of coefficient matrices--ordered consequtively.  Each matrix is
%         expected to be fairly small.
%
%   b   - Array (type `double`) containing the elements of all LLS
%         right-hand sides.
%
%   rsz - Sequence of matrix block row sizes.  The number of matrix blocks
%         contained in 'A' is implicitly assumed to be `numel(rsz)`.
%
%   csz - Sequence of matrix block column sizes.  The number of
%         elements must match `numel(rsz)`.
%
% RETURNS:
%   x - Array (type `double`) containing the elements/entries of the
%       solutions to the sequence of linear least squares problems.
%
% NOTE:
%   The output is algorithmically equivalent to the loop ::
%
%       pA = cumsum([1, rsz .* csz]);
%       pb = cumsum([1, rsz]);
%       px = cumsum([1, csz]);
%       x  = zeros([px(end) - 1, 1]);
%
%       for i = 1 : numel(csz),
%          iA = pA(i) : pA(i + 1) - 1;
%          ib = pb(i) : pb(i + 1) - 1;
%          ix = px(i) : px(i + 1) - 1;
%
%          x(ix) = reshape(A(iA), [rsz(i), csz(i)]) \ b(ib);
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

buildmex('leastsq_svd.c', mwlib('lapack'), mwlib('blas'));

% Now, run it.
[varargout{1:nargout}] = leastsq_svd(varargin{:});
