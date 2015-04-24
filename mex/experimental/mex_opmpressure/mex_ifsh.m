function varargout = mex_ifsh(varargin)
%Discretise and solve flow equation using compiled C code.
%
% SYNOPSIS:
%   [state, wbhp, wflux] = mex_ifsh(state, G, rock, W, bc, src)
%
% DESCRIPTION:
%   Equivalent to the standard MRST call sequence
%
%   fluid = initSingleFluid('mu', 1, 'rho', 1)
%   S     = computeMimeticIP(G, rock)
%   state = solveIncompFlow(state, G, S, fluid, ...
%                           'wells', W, 'src', src, 'bc', bc)
%
%   Note in particular that the inner products are computed at each call.
%
% PARAMETERS:
%   state   - Reservori state.
%
%   G, rock - Grid and rock data structures, respectively.
%
%   W       - Well data structure as defined by 'addWell'.  May be an empty
%             array (interpreted as no wells).
%
%   bc, src - Boundary condition and source data structures as defined by
%             'addBC' and 'addSource', respectively.  Either may be empty.
%
% RETURNS:
%   state   - Updated reservoir state.  Contains new values for
%             'state.pressure' and 'state.flux'.
%
%   wbhp    - Well bottom-hole pressures.  One scalar for each well.  Empty
%             in case of no wells.
%
%   wflux   - Well perforation fluxes.  One scalar value, positive for
%             injection, for each perforation for each flux.  Empty in case
%             of no wells.
%
% SEE ALSO:
%   mex_ip_simple, mex_schur_comp_symm, test_mex_schur_comp_symm.

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

   SRC = { 'mex_ifsh.c', 'mrst_objects.c'             , ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c'), ...
           fullfile(d, '..', 'mrst_api', 'call_umfpack.c') };

   LIBS = { '-lopmpressure', '-lmwumfpack', '-lmwamd', ...
            '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_ifsh(varargin{:});
end
