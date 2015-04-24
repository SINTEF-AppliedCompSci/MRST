function varargout = mex_generate_coarsegrid(varargin)
%Build coarse grid data structure from partition of existing grid (MEX)
%
% SYNOPSIS:
%   CG = mex_generate_coarsegrid(G, p)
%   CG = mex_generate_coarsegrid(G, p, expected_nconn)
%
% PARAMETERS:
%   G - Grid data structure as described in 'grid_structure'.
%
%   p - Partition vector as defined by, e.g., functions 'partitionUI' or
%       'partitionNonUniform'.
%
%   expected_nconn -
%       Number (non-negative integer) of expected fine-scale faces
%       constituting a coarse-scale face.  If expected_nconn==0, then
%       constituent fine-scale faces will not be computed.  On the other
%       hand, if expected_nconn > 0, then constituent fine-scale faces will
%       be derived (similarly to the output of function 'subFaces').  Any
%       positive number may be used, but the implementation is most
%       efficient if 'expected_nconn' is in the same order of magnitude as
%       the typical number of constituent fine-scale faces.
%
%       OPTIONAL.  Default value: expected_nconn=0 (don't compute
%       constituent fine-scale faces (sub-faces)).
%
% RETURNS:
%   CG - Coarse grid data structure as described in 'generateCoarseGrid'.
%        There is, however, a number of subtle differences in the details
%        of this structure as compared to the pure M implementation.
%        Specifically, the coarse faces are numbered differently, and the
%        MEX implementation does not distinguish cardinal directions whence
%        the 'subFaces' function does not produce meaningful results on
%        outer faces.
%
%        If expected_nconn>0, then the 'faces' structure has two additional
%        fields 'subfacePos', and 'subfaces'.  This indirection/data array
%        pair is related such that the sub-faces for coarse face 'i' is
%        located in
%
%            subfaces(subfacePos(i) : subfacePos(i+1) - 1)
%
%        The constituent sub-faces of a particular coarse face may occur in
%        any order.
%
% SEE ALSO:
%   generateCoarseGrid, subFaces.

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

   SRC = { 'mex_generate_coarsegrid.c'                , ...
           fullfile(d, '..', 'mrst_api', 'mrst_api.c')};

   LIBS = { '-lopmpressure', '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_generate_coarsegrid(varargin{:});
end
