function varargout = mex_iteratedJacobiBasisFaster(varargin)
%Undocumented Utility Function

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
