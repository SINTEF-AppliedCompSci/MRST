% INOUT
%   Low-level input routines for ECLIPSE deck.
%
% Files
%   cacheDir          - Define platform-dependent directory for temporary files.
%   convertInputUnits - Convert input data to MRST's strict SI conventions
%   readCache         - Read cached variables associated with given file into caller's workspace.
%   writeCache        - Write callers workspace to matfile in directory ./.cache/
%   writeGRDECL       - Write a GRDECL structure out to permanent file on disk.

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
