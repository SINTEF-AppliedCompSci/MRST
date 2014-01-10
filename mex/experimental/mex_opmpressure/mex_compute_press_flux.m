function varargout = mex_compute_press_flux(varargin)
%Derive pressure and flux from hybrid system using compiled C code.
%
% SYNOPSIS:
%   [v, p] = mex_compute_press_flux(BI, lam, connPos, conns, F, L, q)
%
% PARAMETERS:
%   BI      - Inner product values.  Typically computed using function
%             'mex_ip_simple'.
%
%   lam     - Interface pressure values.  One scalar value for each face in
%             the discretised reservoir model.
%
%   connPos - Indirection map of size [G.cells.num,1] into 'conns' table
%             (i.e., the connections or DOFs).  Specifically, the DOFs
%             connected to cell 'i' are found in the submatrix
%
%                  conns(connPos(i) : connPos(i + 1) - 1)
%
%   conns   - A (connPos(end)-1)-by-1 array of cell connections
%             (local-to-global DOF mapping in FEM parlance).
%
%   F       - Second-to-last return value from 'mex_schur_comp_symm'.
%
%   L       - Last return value from 'mex_schur_comp_symm'.
%
%   q       - Any explicit source terms used to form the original
%             Schur-complement system, or the fifth return value from
%             function 'mex_schur_comp_symm'.
%
% NOTE:
%   The (connPos,conns) array pair is expected to be the output of function
%   'mex_ip_simple'.
%
% RETURNS:
%   v - A SUM(nconn)-by-1 array of half-contact fluxes, ordered by cells.
%
%   p - A NUMEL(nconn)-by-1 array of cell pressure values.
%
% NOTE:
%   This function is the MEX'ed equivalent to the post-processing of
%   function 'schurComplementSymm'.  Note furthermore that this function
%   can only be used in conjunction with function 'mex_schur_comp_symm'.
%
% EXAMPLE:
%   G = computeGeometry(processGRDECL(makeModel3([100, 60, 15])));
%   K = logNormLayers(G.cartDims, [10, 300, 40, 0.1, 100]);
%   rock.perm = bsxfun(@times, [1, 100, 0.1], K(:));
%   rock.perm = convertFrom(rock.perm(G.cells.indexMap, :), ...
%                           milli*darcy);
%
%   [BI, connPos, conns] = mex_ip_simple(G, rock);
%
%   nconn = diff(connPos);
%
%   [i, j] = blockDiagIndex(nconn, nconn);
%   [S, r, F, L] = mex_schur_comp_symm(BI, connPos, conns);
%
%   SS = sparse(double(conns(i)), double(conns(j)), S);
%   R = accumarray(conns, r);
%
%   lam = SS \ R;
%
%   t0 = tic;
%   [v, p] = mex_compute_press_flux(BI, lam, connPos, conns, F, L);
%   toc(t0)
%
% SEE ALSO:
%   mex_ip_simple, mex_schur_comp_symm.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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

   OPTS = { '-g', '-largeArrayDims' };

   SRC = { 'mex_compute_press_flux.c', 'mrst_objects.c', ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c') };

   LIBS = { '-lopmpressure', '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_compute_press_flux(varargin{:});
end
