function varargout = reordersequence_mex(varargin)
% Build MEX edition of same.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   v = version;v = v([1,3]);
   d = fileparts(mfilename('fullpath'));
   CFLAGS = {'CFLAGS="\$CFLAGS', '-O3', '-Wall', '-Wextra', '-ansi', ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align', ...
             '-Wpointer-arith', '-Wbad-function-cast', ...
             '-Wmissing-prototypes ', '-Wstrict-prototypes', ...
             '-Wmissing-declarations', '-Winline', '-Wundef', ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow', ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};

   SRC = {'reordersequence_mex.c', 'reordersequence.c', 'tarjan.c',...
          fullfile(d, '..', '..', '..', 'mrst_api', 'mrst_api.c'),};

   INCLUDE = { ['-I', fullfile(d, '..', '..', '..', 'mrst_api')] };

   OPTS = {'-output', ['reordersequence_mex.', mexext], ...
      '-largeArrayDims', ['-DMATLABVERSION=', v], '-O'};

   buildmex(CFLAGS{:}, INCLUDE{:}, OPTS{:}, SRC{:})

   % Call MEX edition.
   [varargout{1:nargout}] = reordersequence_mex(varargin{:});
end
