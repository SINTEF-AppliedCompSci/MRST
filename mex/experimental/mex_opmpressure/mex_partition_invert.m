function varargout = mex_partition_invert(varargin)
%Invert cell-to-block map (creating block-to-cell) using compiled C code.
%
% SYNOPSIS:
%   [pb2c, b2c]      = mex_partition_invert(p)
%   [pb2c, b2c, loc] = mex_partition_invert(p)
%
% PARAMETERS:
%   p - Partition vector.  Should not contain any empty blocks.  Use
%       function 'mex_partition_compress' to remove empty blocks/bins.
%
% RETURNS:
%   pb2c - Indirection map of size [MAX(p) + 1, 1] into the 'b2c' map
%          array.  Specifically, the cells of block 'b' are stored in
%
%                b2c(pb2c(b) : pb2c(b + 1) - 1)
%
%   b2c  - Inverse cell map.  The entries in pb2c(b):pb2c(b+1)-1 correspond
%          to the result of FIND(p == b).
%
%   loc  - Local index within a block/bin.  Specifically,
%
%             loc(i) == FIND(b2c == i) - pb2c(p(i)) + 1
%
%          OPTIONAL.  Only returned (and computed) if specifically
%          requested.
%
% SEE ALSO:
%   mex_partition_ui, mex_partition_compress.

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

   OPTS = { '-O', '-largeArrayDims' };

   SRC = { 'mex_partition_invert.c' };

   LIBS = { '-lopmpressure', '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_partition_invert(varargin{:});
end
