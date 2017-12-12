function [G,F,fracture] = gridFracture2D(G,fracture,varargin)
% gridFracture2D imposes a fracture grid given the underlying matrix grid,
% information about the fracture lines and desired grid resolution
% (optional).
%
% SYNOPSIS:
%   [G,F,fracture] = gridFracture2D(G, fracture)
%   [G,F,fracture] = gridFracture2D(G, fracture, 'pn1', pv1)
%
% REQUIRED PARAMETERS:
%
%   G         - Matrix grid structure with the sub-structure
%               G.cells.fracture. See processFracture2D.
%
%   fracture  - Structure containing information about fracture networks,
%               independent fractures and their conductivity towards the
%               matrix cells they penetrate. See processFracture2D.
%
% OPTIONAL PARAMETERS:
%   Same as assembleFracNodes2D
%
% RETURNS:
%   G        - Matrix grid structure with sub-structure G.FracGrid (see
%              FracTensorGrid2D).
%
%   F        - Same as assembleFracNodes2D.
%
%   fracture - structure with the added structure - 'intersections' which
%              contains the following fields:
%              (a) lines - n-by-2 matrix of pairs of intersecting lines
%              where n is the total number of fracture intersections.
%              (b) coords - coordinates of intersection correspondng to
%              field lines
%
% NOTE: 
%   This function calls assembleFracNodes2D and FracTensorGrid2D
%   internally
%
% SEE ALSO:
%   assembleFracNodes2D, FracTensorGrid2D, processFracture2D,
%   getIndepNetwork, markcells2D, tensorGrid

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


opt = struct('verbose'     , mrstVerbose,...
             'assemblyType', 1,...
             'min_size'    , 0.1,...
             'cell_size'   , 1);
opt = merge_options(opt, varargin{:});


% Divide fracture lines
dispif(opt.verbose, 'Gridding fracture lines...\n\n');
[F,fracture] = assembleFracNodes2D(G,fracture,...
               'assemblyType',opt.assemblyType,...
               'min_size',opt.min_size,'cell_size',opt.cell_size);


% Build fracture grid
dispif(opt.verbose, 'Building fracture grid as tensorGrid...\n\n');
G = FracTensorGrid2D(G,F,fracture.aperture);

return