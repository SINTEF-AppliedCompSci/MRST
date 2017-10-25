function varargout = amgcl_matlab(varargin)
% MEX-gateway to AMGCL solver.
% For more information about AMGCL, please see:
% Documentation: http://amgcl.readthedocs.io/en/latest/
% GitHub repository: https://github.com/ddemidov/amgcl/tree/master/examples
%
% Tested with commit c1c1fc565f96c1bb3adc0a1bf92429ad744036b4.

   global AMGCLPATH
   global BOOSTPATH
   
   CFLAGS = {'CXXFLAGS="$CXXFLAGS', '-std=c++11', '-O3', '-fopenmp"'};
   LDFLAGS = {'LDFLAGS="$LDFLAGS', '-fopenmp"'};

   INCLUDE = {['-I', BOOSTPATH], ['-I', AMGCLPATH]};

   OPTS = {};

   SRC = {'amgcl_matlab.cpp'};

   [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params;

   
   buildmex(CFLAGS{:}, CXXFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});
        
   % Call MEX'ed edition.
   [varargout{1:nargout}] = amgcl_matlab(varargin{:});
end

function [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params
   e = mexext('all');
   a = e(strcmp({ e.ext }, mexext)).arch;

   if ispc

      mwlib = @(lib) ...
      fullfile(matlabroot, 'extern', 'lib', a, 'microsoft', ...
               ['libmw', lib, '.lib']);

      CXXFLAGS  = { };
      LINK      = { };
      iomp5     = { };
      libstdcpp = { };

   elseif isunix

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
