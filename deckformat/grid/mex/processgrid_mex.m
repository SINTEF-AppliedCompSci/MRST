function varargout = processgrid_mex(varargin)
%Compute grid topology and geometry from pillar grid description.
%
% SYNOPSIS:
%   G = processGRDECL(grdecl)
%
% PARAMETERS:
%   grdecl - Raw pillar grid structure, as defined by function
%            'readGRDECL', with fields COORDS, ZCORN and, possibly, ACTNUM.
%
% RETURNS:
%   G      - Valid grid definition containing connectivity, cell
%            geometry, face geometry and unique nodes.
%
% EXAMPLE:
%   G = processgrid(readGRDECL('small.grdecl'));
%   plotGrid(G); view(10,45);
%
% SEE ALSO:
%   `processGRDECL`, `readGRDECL`, `deactivateZeroPoro`, `removeCells`.

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


% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% Mex gateway by Jostein R. Natvig, SINTEF ICT.


   % Build MEX edition of same.
   %

   v = version;v = v([1,3]);

   CFLAGS = {'CFLAGS="$CFLAGS', '-O3', '-Wall', '-Wextra', '-ansi', ...
             '-pedantic', '-Wformat-nonliteral',  '-Wcast-align', ...
             '-Wpointer-arith', '-Wbad-function-cast', ...
             '-Wmissing-prototypes', '-Wstrict-prototypes', ...
             '-Wmissing-declarations', '-Winline', '-Wundef', ...
             '-Wnested-externs', '-Wcast-qual', '-Wshadow', ...
             '-Wconversion', '-Wwrite-strings', '-Wno-conversion', ...
             '-Wchar-subscripts', '-Wredundant-decls"'};

   SRC = {'processgrid.c', 'preprocess.c', 'uniquepoints.c', ...
          'facetopology.c', 'mxgrdecl.c'};

   INCLUDE = {};

   OPTS = {'-output', 'processgrid_mex', ...
           '-largeArrayDims', ['-DMATLABVERSION=', v], '-O'};

   buildmex(CFLAGS{:}, INCLUDE{:}, OPTS{:}, SRC{:})

   % Call MEX edition.
   [varargout{1:nargout}] = processgrid_mex(varargin{:});
end

