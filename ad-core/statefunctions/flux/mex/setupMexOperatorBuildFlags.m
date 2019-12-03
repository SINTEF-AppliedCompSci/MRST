function [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags()
%Undocumented Utility Function

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

   a = computer('arch');

   if ispc
      mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, ...
               'microsoft', ['libmw', lib, '.lib']);

      % Note explicit /EHsc to enable C++ exception handling
      CXXFLAGS  = { 'COMPFLAGS=/EHsc /MD /openmp /O2' };
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
