% Supporting utilities for MATLAB Reservoir Simulation Toolbox (MRST).
%
% Files
%   ADI.m                   - ADI class: simple implementation of automatic differentiation for easy construction of jacobian matrices.
%   ROOTDIR.m               - Retrieve full path of Toolbox installation directory.
%   blockDiagIndex.m        - Compute subscript or linear index to nonzeros of block-diagonal matrix
%   buildmex.m              - Wrapper around MEX which abstracts away details of pathname generation.
%   cellFlux2faceFlux.m     - Transform cell-based flux field to face-based.
%   dinterpq1.m             - Compute derivative of piecewise linear interpolant.
%   dispif.m                - Produce textual output contingent upon predicate.
%   faceFlux2cellFlux.m     - Transform face-based flux field to cell-based.
%   faceFlux2cellVelocity.m - Transform face-based flux field to one constant velocity per cell.
%   findFilesSubfolders.m   - Search for files in a directory hierarchy
%   formatTimeRange.m       - Small utility which returns a human readable string from seconds.
%   geomspace.m             - Geometrically spaced vector.
%   initVariablesADI.m      - Initialize a set of automatic differentiation variables
%   invv.m                  - INVV(A, sz)
%   mcolon.m                - Compute vector of consecutive indices from separate lower/upper bounds.
%   md5sum.m                - md5sum - Compute md5 check sum of all input arguments
%   md5sum_fallback.m       - Alternative implementation of md5sum for systems without C compiler
%   merge_options.m         - Override default control options.
%   mrstDebug.m             - Globally control default settings for MRST debugging information.
%   mrstModule.m            - Query or modify list of activated add-on MRST modules
%   mrstVerbose.m           - Globally control default settings for MRST verbose information.
%   msgid.m                 - Construct Error/Warning message ID by prepending function name.
%   require.m               - Announce and enforce module dependency.
%   rldecode.m              - Decompress run length encoding of array A along dimension dim.
%   rlencode.m              - Compute run length encoding of array A along dimension dim.
%   ticif.m                 - Evaluate function TIC if input is true.
%   tocif.m                 - Evaluate function TOC if input is true.

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
