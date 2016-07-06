function varargout = mex_partition_process(varargin)
%Split disconnected blocks into new blocks using compiled C code.
%
% SYNOPSIS:
%   p = mex_partition_process(p, Neighbours)
%
% PARAMETERS:
%   p - Partition vector.  Should not contain any empty blocks.  Use
%       function 'mex_partition_compress' to remove empty blocks/bins.
%
%   Neighbours -
%       Neighbour definition.  An m-by-two numeric array of cell-to-cell
%       connections.  Often equal to the 'faces.neighbors' array of an MRST
%       'grid_structure', but general definitions are supported.
%
% RETURNS:
%   p - Updated partition vector.  Disconnected blocks of original partition
%       vector are assigned new block numbers.  Specifically, if block 'b'
%       contains 'ncomp' separate connected components, then ncomp-1 new
%       blocks are created from block 'b'.
%
% NOTE:
%   This implementation corresponds to setting Reconnect=FALSE in the
%   'processPartition' function when 'Neighbour' is G.faces.neighbor of a
%   grid_structure, 'G'.
%
% SEE ALSO:
%   processPartition.

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

   SRC = { 'mex_partition_process.c' };

   LIBS = { '-lopmpressure', '-lmwlapack', '-lmwblas' };

   buildmex(CFLAGS{:}, LDFLAGS{:}, ...
            INCLUDE{:}, LINK{:}, OPTS{:}, SRC{:}, LIBS{:});

   % Call MEX'ed edition.
   [varargout{1:nargout}] = mex_partition_process(varargin{:});
end
