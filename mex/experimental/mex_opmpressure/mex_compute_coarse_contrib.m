function varargout = mex_compute_coarse_contrib(varargin)
%Generate contributions to coarse-scale linsys using compiled C code.
%
% SYNOPSIS:
%   [cell_ip, Binv] = ...
%        mex_compute_coarse_contrib(BIf, Psi, p, pconn, dof_pos, blk_ncf)
%
% PARAMETERS:
%   BIf     - Fine-scale (inverse) inner product matrices.
%
%   Psi     - Basis function values.  Non-zeros only.  Ordered by blocks.
%
%   p       - Partition vector.
%
%   pconn   - Fine-scale indirection array into connection table.
%             Typically corresponds to G.cells.facePos.
%
%   dof_pos - Coarse-scale indirection array into (coarse) connection
%             table.  Roughly equivalent to CG.cells.facePos, but only
%             for active faces.
%
%   blk_ncf - Number of (fine-scale) half-contacts per block.
%
% RETURNS:
%   cell_ip - Cell contributions to coarse-scale (block) inner products.
%
%   Binv    - Coarse-scale inverse inner product.  Semantically
%             equivalent to 'BI' output of 'mex_ip_simple', but for the
%             coarse grid represented by 'p'.
%
% NOTE:
%   This function is only intended for testing/developing a compiled
%   language implementation of the MsMFE method.
%
% SEE ALSO:
%   mex_ip_simple.

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

   SRC = { 'mex_compute_coarse_contrib.c' };

   LIBS = { '-lopmpressure', '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_generate_coarse_contrib(varargin{:});
end
