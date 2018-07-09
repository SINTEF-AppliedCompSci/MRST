function varargout = amgcl_matlab(varargin)
%MEX-gateway to AMGCL Linear Solver Software
%
% For more information about AMGCL, please see:
% Documentation: http://amgcl.readthedocs.io/en/latest/
% GitHub repository: https://github.com/ddemidov/amgcl/tree/master/examples
%
% SYNOPSIS:
%    x              = amgcl_matlab(A, b, amg_opt, tol, maxIter, id)
%   [x, err]        = amgcl_matlab(...)
%   [x, err, nIter] = amgcl_matlab(...)
%
% PARAMETERS:
%   A       - Sparse coefficient matrix of system of simultaneous linear
%             equations.
%
%   b       - System right-hand side.
%
%   amg_opt - AMGCL options structure as defined by function
%             `getAMGCLMexStruct`.
%
%   tol     - Relative residual reduction tolerance.  Positive scalar.
%
%   maxIter - Maximum number of linear iterations.
%
%   id      - Solver method ID.  Integer.  Supported values are `1` for the
%             regular solver and `2` for the CPR solver.
%
% RETURNS:
%   x     - Solution vector.
%
%   err   - Norm of residual at end of solution process.
%
%   nIter - Number of linear iterations.
%
% NOTE:
%   This gateway was last tested with commit
%
%      946398e535cf2586ca59f37eaf8aa9e72f43fde4
%
%   of the AMGCL software (from GitHub).
%
%   * Linux (Mint 17.2):
%       GCC   4.9.4
%       Boost 1.63.0
%
%   * MS Windows 10 (1607):
%       MSVC  19.14.26430 (Visual Studio 15.7)
%       Boost 1.67.0
%
% SEE ALSO:
%   `callAMGCL`, `getAMGCLMexStruct`.

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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

   if ~exist('verLessThan', 'file') || verLessThan('matlab', '8.3.0')
      error(['Automated Build Script for ''amgcl_matlab'' is not ', ...
             'Supported in MATLABs prior to 8.3.0 (R2014a)']);
   end

   global AMGCLPATH
   global BOOSTPATH

   if ~valid_global_path(AMGCLPATH)
      error(['Cannot Build AMGCL MEX Gateway Unless GLOBAL ', ...
             '''AMGCLPATH'' Variable is Set in Current MATLAB Session']);
   end

   if ~valid_global_path(BOOSTPATH)
      error(['Cannot Build AMGCL MEX Gateway Unless GLOBAL ', ...
             '''BOOSTPATH'' Variable is Set in Current MATLAB Session']);
   end

   INCLUDE = strcat('-I', { BOOSTPATH, AMGCLPATH });

   OPTS = { '-O' };

   SRC = {'amgcl_matlab.cpp'};

   [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = amgcl_matlab(varargin{:});
end

%--------------------------------------------------------------------------

function [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params
   a = computer('arch');

   if ispc

      mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, ...
               'microsoft', ['libmw', lib, '.lib']);

      % Note explicit /EHsc to enable C++ exception handling
      CXXFLAGS  = { 'COMPFLAGS=/EHsc /MD /openmp' };
      LINK      = { ['-L', fullfile(matlabroot, 'bin', a) ]};
      iomp5     = { ['LINKFLAGS=$LINKFLAGS ', ...
                     '/nodefaultlib:vcomp libiomp5md.lib' ]};
      libstdcpp = {};

   elseif isunix

       mwlib = @(lib) ['-lmw', lib];

       CXXFLAGS = ...
          { ['CXXFLAGS=$CXXFLAGS -D_GNU_SOURCE -fPIC -O3 ', ...
             '-std=c++11 -ffast-math -march=native -fopenmp'] };

       LINK = { ['-L', fullfile(matlabroot, 'sys', 'os', a)] };

       iomp5     = { '-liomp5', 'LDFLAGS=$LDFLAGS -fopenmp' };
       libstdcpp = { '-lstdc++' };
   end

   LIBS = [ iomp5, { mwlib('lapack'), mwlib('blas') }, libstdcpp ];
end

%--------------------------------------------------------------------------

function tf = valid_global_path(p)
   tf = ~isempty(p) && ischar(p) && isdir(p);
end
