function [CXXFLAGS, LINK, LIBS] = setupMexOperatorBuildFlags()
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
