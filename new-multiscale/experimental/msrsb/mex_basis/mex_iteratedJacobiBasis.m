function varargout = mex_iteratedJacobiBasis(varargin)
   d = fileparts(mfilename('fullpath'));

   CFLAGS = {'CFLAGS="\$CFLAGS', '-fopenmp', '-Wall','-O2', '-Wextra', '-ansi'     , ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align' , ...
             '-Wpointer-arith', '-Wbad-function-cast'            , ...
             '-Wmissing-prototypes', '-Wstrict-prototypes'      , ...
             '-Wmissing-declarations', '-Winline', '-Wundef'     , ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow'       , ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-fopenmp', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};

   LDFLAGS = {'-lgomp'};

   INCLUDE = { ['-I', fullfile(d, '..', 'mrst_api')] };

   LINK = {  };

   OPTS = { '-O', '-largeArrayDims'};

   SRC = {'mex_iteratedJacobiBasis.c', 'jacobi_basis.c'};

   LIBS = {  };
    
   
   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_iteratedJacobiBasis(varargin{:});
end