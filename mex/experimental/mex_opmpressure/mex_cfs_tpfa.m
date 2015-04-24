function varargout = mex_cfs_tpfa(varargin)
%Perform single fixed point iteration using compiled C code.
%
% SYNOPSIS:
%   [p, Tf, Gf] = mex_cfs_tpfa(G, T, fq, Zf, dt, p0, pvol, bc, src)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   T       - Face (not h-face) transmissibility.  Unaffected by mobility.
%
%   fq      - Fluid quantitity structure.  Must contain the following
%             fields:
%                - np    - Number of phases (scalar)
%                - ct    - Total compressibility (per cell)
%                - Lt    - Total mobility (per cell)
%                - A.c   - (R / B) per cell (np-by-np*nc num. array)
%                - A.f   - (R / B) per face (np-by-np*nf num. array)
%                - pmobf - Phase mobility per face in an np-by-nf
%                          numeric array
%
%   Zf      - Gravity/capillary pressure contributions.  np-by-nf array
%
%   dt      - Time step.  Positive scalar.
%
%   p0      - Pressure field at previous pressure step.  One positive
%             scalar for each cell.
%
%   pvol    - Pore volume.  One positive scalar for each cell.
%
%   bc, src - Boundary condition and source data structures as defined by
%             'addBC' and 'addSource', respectively.  Either may be empty.
%
% RETURNS:
%   p  - New pressure values.
%
%   Tf - Face transmissibility.  T .* (R/B)_c\(R/B)_f * [pmobf].
%
%   Gf - Gravity flux contributions.  T .* (R/B)_c\(R/B)_f * [Zf].
%
% SEE ALSO:
%   matrixBlocksFromSparse.

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

   if ismac,
      LDFLAGS = {'LDFLAGS="\$LDFLAGS', ...
                 '-Wl,-rpath', ['-Wl,', fullfile(d, '..', 'lib'), '"']};
   else
      assert (isunix, ['Function is not supported on platforms other ' ...
                       'than MacOS or Unix']);
      LDFLAGS = {'LDFLAGS="\$LDFLAGS', ...
                 ['-Wl,-rpath=', fullfile(d, '..', 'lib'), '"']};
   end

   INCLUDE = { '-I/usr/include/suitesparse'                          , ...
               '-I/usr/local/include/SuiteSparse'                    , ...
              ['-I', fullfile(d, '..', 'include', 'opm', 'pressure')], ...
              ['-I', fullfile(d, '..', 'mrst_api')] };

   LINK = { ['-L', fullfile(d, '..', 'lib')] };

   OPTS = { '-O', '-largeArrayDims' };

   SRC = { 'mex_cfs_tpfa.c', 'mrst_objects.c'         , ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c'), ...
           fullfile(d, '..', 'mrst_api', 'call_umfpack.c') };

   LIBS = { '-lopmpressure', '-lmwumfpack', '-lmwamd', ...
            '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_cfs_tpfa(varargin{:});
end
