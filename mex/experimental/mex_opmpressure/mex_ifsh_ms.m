function varargout = mex_ifsh_ms(varargin)
%Discretise and solve flow equation using compiled C code.
%
% SYNOPSIS:
%   state = mex_ifsh_ms(state, G, rock, p, src)
%
% PARAMETERS:
%   state   - Reservoir state.
%
%   G, rock - Grid and rock data structures, respectively.
%
%   p       - Partition vector.
%
%   src     - Source data structures as defined by 'addSource'.
%
% RETURNS:
%   state   - Updated reservoir state.  Contains new values for
%             'state.pressure' and 'state.flux'.
%
% SEE ALSO:
%   mex_ifsh, mex_generate_coarsegrid, mex_generate_coarsesystem.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   d = fileparts(mfilename('fullpath'));

   CFLAGS = {'CFLAGS="\$CFLAGS', '-Wall', '-Wextra', '-ansi', ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align', ...
             '-Wpointer-arith', '-Wbad-function-cast', ...
             '-Wmissing-prototypes', '-Wstrict-prototypes', ...
             '-Wmissing-declarations', '-Winline', '-Wundef', ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow', ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};

   LDFLAGS = {'LDFLAGS="\$LDFLAGS', ...
              ['-Wl,-rpath=', fullfile(d, '..', 'lib'), '"']};

   INCLUDE = { '-I/usr/include/suitesparse'                          , ...
              ['-I', fullfile(d, '..', 'include', 'opm', 'pressure')], ...
              ['-I', fullfile(d, '..', 'mrst_api')] };

   LINK = { ['-L', fullfile(d, '..', 'lib')] };

   OPTS = { '-O', '-largeArrayDims' };

   SRC = { 'mex_ifsh_ms.c'                            , ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c'), ...
           fullfile(d, '..', 'mrst_api', 'call_umfpack.c') };

   LIBS = { '-lopmpressure', '-lmwumfpack', '-lmwamd', ...
            '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_ifsh_ms(varargin{:});
end
