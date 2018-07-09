function varargout = mex_cppMultiscaleBasis(varargin)
% Internal build routine for MsRSB accelerated basis functions

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

   CFLAGS = {'CXXFLAGS="$CXXFLAGS', '-std=c++11', '-O3', '-fopenmp"'};
   LDFLAGS = {'LDFLAGS="$LDFLAGS', '-fopenmp"'};

   INCLUDE = {};

   OPTS = { '-O', '-largeArrayDims', '-DUSEMEX=""'};

   SRC = {'main.cpp', 'basis_solver.cpp'};

   [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params;

   
   buildmex(CFLAGS{:}, CXXFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});
        
   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_cppMultiscaleBasis(varargin{:});
end

function [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params
   e = mexext('all');
   a = e(strcmp({ e.ext }, mexext)).arch;

   if ispc,

      mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, 'microsoft', ...
               ['libmw', lib, '.lib']);

      CXXFLAGS  = { };
      LINK      = { };
      iomp5     = { };
      libstdcpp = { };

   elseif isunix,

       mwlib = @(lib) ['-lmw', lib];

       CXXFLAGS = ...
          {'CXXFLAGS="-D_GNU_SOURCE', '-fPIC', '-O3', '-std=c++0x',  ...
           '-fopenmp"'};

       LINK = { ['-L', fullfile(matlabroot, 'sys', 'os', a)] };

       iomp5     = { '-liomp5' };
       libstdcpp = { '-lstdc++' };
   end

   LIBS = [ iomp5, { mwlib('lapack'), mwlib('blas') }, libstdcpp ];
end
