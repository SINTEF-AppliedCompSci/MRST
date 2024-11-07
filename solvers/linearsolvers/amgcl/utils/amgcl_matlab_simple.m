function varargout = amgcl_matlab_simple(varargin)
%MEX-gateway to AMGCL Linear Solver Software
%
% For more information about AMGCL, please see:
% Documentation: http://amgcl.readthedocs.io/en/latest/
% GitHub repository: https://github.com/ddemidov/amgcl/tree/master/examples
%
% SYNOPSIS:
%    x              = amgcl_matlab(A, b, amg_opt, id)
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
%   For first-time use, this gateway will attempt to build a MEX-file using
%   the configured C++ compiler. In order to do this, the paths to the
%   dependencies AMGCL and BOOST must be specified via global variables.
%
%   For instance, this can be done by executing the following:
%   global BOOSTPATH AMGCLPATH
%   % Define paths
%   BOOSTPATH = '/path/to/boost';
%   AMGCLPATH = '/path/to/amgcl-repo';
%   % Build the binaries
%   amgcl_matlab_simple();
%
%   In addition, you will require a working C++ compiler (see mex -setup).
%
%   This gateway was last tested with commit
%
%      91348b025144b61a2d3b7417988373d0dc8c8d00
%
%   of the AMGCL software (from GitHub: https://github.com/ddemidov/amgcl).
%
%   * Linux (Mint 17.2):
%       GCC   4.9.4
%       Boost 1.63.0
%
%   * MS Windows 10 (1607):
%       MSVC  19.16.27024.1 (Visual Studio 15.9)
%       Boost 1.68.0
%
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

   [AMGCLPATH, BOOSTPATH] = getAMGCLDependencyPaths();
   if ~valid_global_path(AMGCLPATH)
      error(['Cannot Build AMGCL MEX Gateway Unless GLOBAL ', ...
             '''AMGCLPATH'' Variable is Set in Current MATLAB Session.', ...
             ' For detailed build instructions, see "help amgcl_matlab".']);
   end

   if ~valid_global_path(BOOSTPATH)
      error(['Cannot Build AMGCL MEX Gateway Unless GLOBAL ', ...
             '''BOOSTPATH'' Variable is Set in Current MATLAB Session.', ...
             ' For detailed build instructions, see "help amgcl_matlab".']);
   end

   INCLUDE = strcat('-I', { BOOSTPATH, AMGCLPATH });

   OPTS = { '-O' };

   SRC = {'amgcl_matlab_simple.cpp'};

   [CXXFLAGS, LINK, LIBS] =  mrstDefaultMexFlags('AMGCL_ASYNCSETUP');

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = amgcl_matlab_simple(varargin{:});
end

%--------------------------------------------------------------------------

function tf = valid_global_path(p)
   tf = ~isempty(p) && ischar(p) && isdir(p);
end
