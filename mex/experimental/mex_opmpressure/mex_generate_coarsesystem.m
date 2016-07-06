function varargout = mex_generate_coarsesystem(varargin)
%Build coarse system data structure from partition of existing grid (MEX)
%
% SYNOPSIS:
%   CS = mex_generate_coarsesystem(G, rock, CG, Lt, src)
%
% PARAMETERS:
%   G    - Grid data structure as described in 'grid_structure'.
%
%   rock - Rock data structure with valid field 'perm'.
%
%   CG   - Coarse-grid stucture as defined by function
%          'mex_generate_coarsegrid'.  Must, in particular, support the
%          fields 'faces.subfacePos' and 'faces.subfaces'
%
%   Lt   - Vector of total mobilities.  One (positive) scalar for each cell
%          in the underlying fine-scale grid, 'G'.
%
%   src  - Source data structure as defined by function 'addSource'.  Pass
%          an empty array if there are no explicit source terms in the
%          problem.
%
% RETURNS:
%   CS - Coarse linear system data structure similar to that returned by
%        'generateCoarseSystem'.  The structure has the following fields:
%
%            - Dof2Conn  -- Map DOFs to coarse connections
%            - blkDofPos -- Start pointers to each block's DOFs
%            - basisPos  -- Start pointers to each block's BFs
%            - cellIPPos -- Start pointers to each block's cell IPs
%            - blkDof    -- Each block's DOFs
%            - basis     -- Each block's BFs
%            - cellIP    -- Each block's fine-scale cell IPs
%            - BI        -- Each block's coarse-scale inverse B (unset)
%
%        Note: The start pointer and block-DOF arrays are zero-based.  Use
%        [ I(x)+1 : I(x+1) ] to index into the corresponding data array.
%
% SEE ALSO:
%   mex_generate_coarsegrid, subFaces.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

   SRC = { 'mex_generate_coarsesystem.c', 'mrst_objects.c' , ...
           fullfile(d, '..', 'mrst_api', 'call_umfpack.c') , ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c')     };

   LIBS = { '-lopmpressure', '-lmwumfpack', '-lmwamd', ...
            '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_generate_coarsesystem(varargin{:});
end
