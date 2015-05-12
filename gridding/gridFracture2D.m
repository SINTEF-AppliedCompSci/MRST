function [G,F,fracture] = gridFracture2D(G,fracture,varargin)
% gridFracture2D grids a fracture given the matrix grid and fracture lines.
%
% SYNOPSIS:
%   [G,F,fracture] = gridFracture2D(G, fracture)
%   [G,F,fracture] = gridFracture2D(G, fracture, 'pn1', pv1)
%
% REQUIRED PARAMETERS:
%
%   G         - Matrix grid structure with the sub-structure
%               G.cells.fracture (see markcells).
%
%   fracture  - Structure containing information pertaining to independant
%               fracture networks, individual fracture lines and matrix
%               cells with embedded fractures. See processFracture.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   Same as assembleFracNodes2D
%
% RETURNS:
%   G  - Matrix grid structure with structure G.FracGrid (see
%        FracTensorGrid2D)
%
%   F  - Same as assembleFracNodes2D
%
% NOTE: 
%   This function calls assembleFracNodes2D and FracTensorGrid2D
%   internally
%
% SEE ALSO:
%   assembleFracNodes2D, FracTensorGrid2D, processFracture,
%   getIndepNetwork, markcells 

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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


opt = struct('verbose', false,'assemblyType',1,'min_size',0.1,'cell_size',1);
opt = merge_options(opt, varargin{:});


% Divide fracture lines
dispif(opt.verbose, 'Gridding fracture lines...\n\n');
[F,fracture] = assembleFracNodes2D(G,fracture,'assemblyType',opt.assemblyType,...
    'min_size',opt.min_size,'cell_size',opt.cell_size);


% Build fracture grid
dispif(opt.verbose, 'Building fracture grid as tensorGrid...\n\n');
G = FracTensorGrid2D(G,F,fracture.aperture);

return