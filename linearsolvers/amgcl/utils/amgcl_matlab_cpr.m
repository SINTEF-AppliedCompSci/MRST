function varargout = amgcl_matlab_cpr(varargin)
% MEX-gateway to AMGCL solver.
% For more information about AMGCL, please see:
% Documentation: http://amgcl.readthedocs.io/en/latest/
% GitHub repository: https://github.com/ddemidov/amgcl/tree/master/examples
%
% Tested with commit c1c1fc565f96c1bb3adc0a1bf92429ad744036b4.

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

   SRC = {'amgcl_matlab_cpr.cpp'};

   [CXXFLAGS, LINK, LIBS] = setup_machdep_build_params();

   buildmex(OPTS{:}, INCLUDE{:}, CXXFLAGS{:}, SRC{:}, LINK{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = amgcl_matlab_cpr(varargin{:});
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
