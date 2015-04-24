function varargout = mex_schur_comp(varargin)
%Compute hybrid system component matrices using compiled C code.
%
% SYNOPSIS:
%   [S, r, F1, F2, L] = mex_schur_comp(BI, connPos, conns)
%
% PARAMETERS:
%   BI      - Inner product values.
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
% NOTE:
%   The (connPos,conns) array pair is expected to be the output of function
%   'mex_ip_simple'.
%
% RETURNS:
%   S  - A SUM(nconn .^ 2)-by-1 array of unassembled system matrix values,
%        ordered by cells.
%
%   r  - A SUM(nconn)-by-1 array of unassemble system rhs values, ordered
%        by cells.
%
%   F1 - A SUM(nconn)-by-1 array of C'*inv(B) values, ordered by cells.
%
%   F2 - A SUM(nconn)-by-1 array of C'*inv(B) values, ordered by cells.
%
%   L  - A G.cells.num-by-1 array of C'*inv(B)*C values, ordered by cells.
%
% SEE ALSO:
%   mex_ip_simple, mex_compute_press_flux.

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

   SRC = { 'mex_schur_comp.c', ...
           fullfile(d, '..', 'mrst_api', 'call_umfpack.c') };

   LIBS = { '-lopmpressure', '-lmwumfpack', '-lmwamd', ...
            '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_schur_comp_symm(varargin{:});
end
