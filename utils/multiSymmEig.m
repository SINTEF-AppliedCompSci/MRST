function varargout = multiSymmEig(varargin)
%Solve sequence of symmetric eigenvalue problems using LAPACK
%
% SYNOPSIS:
%    d     = multiSymmEig(A, sz)
%   [d, v] = multiSymmEig(A, sz)
%
% PARAMETERS:
%   A  - Array (type `double`) containing the elements/entries of a sequence
%        of coefficient matrices--ordered consequtively.  Each matrix is
%        expected to be square, symmetric and fairly small.
%
%   sz - Sequence of matrix block sizes.  The number of matrix blocks
%        contained in `A` is implicitly assumed to be `numel(rsz)`.
%
% RETURNS:
%   d - Array (type `double`) containing the eigenvalues of the individual
%       eigenvalue problems.
%
%   v - Array (type `double`) of eigenvectors.  Optional.  Only returned if
%       specifically requested.
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

   [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params;

   INCLUDE = { };

   OPTS = { '-O', '-largeArrayDims' };

   SRC = { 'multiSymmEig.cpp' };

   % Build MEX file
   buildmex(CXXFLAGS{:}, INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Now, run it.
   [varargout{1:nargout}] = multiSymmEig(varargin{:});
end

%--------------------------------------------------------------------------

function [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params
   a = computer('arch');

   if ispc,

      mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, 'microsoft', ...
               ['libmw', lib, '.lib']);

      CXXFLAGS  = { '/MD', '/openmp' };
      LINK      = { };
      iomp5     = { '/link', '/nodefaultlib:vcomp', 'libiomp5md' };
      libstdcpp = { };

   elseif isunix,

       mwlib = @(lib) ['-lmw', lib];

       CXXFLAGS = ...
          {'CXXFLAGS="-D_GNU_SOURCE', '-fPIC', '-O3', '-std=c++0x',  ...
           '-fopenmp', '-Wall', '-Wextra', '-pedantic',              ...
           '-Wformat-nonliteral', '-Wcast-align', '-Wpointer-arith', ...
           '-Wmissing-declarations', '-Wundef', '-Wcast-qual',       ...
           '-Wshadow', '-Wwrite-strings', '-Wchar-subscripts',       ...
           '-Wredundant-decls"'};

       LINK = { ['-L', fullfile(matlabroot, 'sys', 'os', a)] };

       iomp5     = { '-liomp5' };
       libstdcpp = { '-lstdc++' };
   end

   LIBS = [ iomp5, { mwlib('lapack'), mwlib('blas') }, libstdcpp ];
end
