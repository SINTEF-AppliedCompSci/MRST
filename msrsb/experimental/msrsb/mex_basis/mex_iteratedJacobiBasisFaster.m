function varargout = mex_iteratedJacobiBasisFaster(varargin)
   d = fileparts(mfilename('fullpath'));
    
   if 1
   CFLAGS = {'CFLAGS="\$CFLAGS', '-fopenmp', '-Wall','-O3', '-Wextra', '-ansi'     , ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align' , ...
             '-Wpointer-arith', '-Wbad-function-cast'            , ...
             '-Wmissing-prototypes', '-Wstrict-prototypes'      , ...
             '-Wmissing-declarations', '-Winline', '-Wundef'     , ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow'       , ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-fopenmp', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};
   else
   CFLAGS = {'CFLAGS="\$CFLAGS', '-fopenmp"'};
   end
   LDFLAGS = {'-lgomp'};

   INCLUDE = { ['-I', fullfile(d, '..', 'mrst_api')] };

   LINK = {  };

   OPTS = { '-O', '-largeArrayDims'};

   SRC = {'mex_iteratedJacobiBasisFaster.c', 'jacobi_basis_faster.c'};

   LIBS = {  };
    
   
   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_iteratedJacobiBasisFaster(varargin{:});
end