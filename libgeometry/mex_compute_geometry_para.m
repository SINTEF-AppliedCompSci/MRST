function varargout = mex_compute_geometry(varargin)
   d = fileparts(mfilename('fullpath'));

   CFLAGS = {'CFLAGS="\$CFLAGS','-O3', '-fopenmp', '-Wall', '-Wextra', '-ansi'     , ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align' , ...
             '-Wpointer-arith', '-Wbad-function-cast'            , ...
             '-Wmissing-prototypes', '-Wstrict-prototypes'      , ...
             '-Wmissing-declarations', '-Winline', '-Wundef'     , ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow'       , ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};

   LDFLAGS = {'LDFLAGS="\$LDFLAGS','-fopenmp"'};

   INCLUDE = { ['-I', fullfile(d, '..', 'mrst_api')] };

   LINK = { };

   OPTS = { '-O', '-largeArrayDims'};

   SRC = { 'mex_compute_geometry_para.c', 'geometry_para.c', ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c') };

   LIBS = {};

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_geometry_para(varargin{:});
end
