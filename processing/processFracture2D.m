function [G,fracture] = processFracture2D(G,fl,varargin)
% processFracture extracts independant fracture networks and stores matrix
% cells containing fractures given a matrix grid and a set of fracture
% lines
%
% % SYNOPSIS:
%   [G,F,fracture] = gridFracture2D(G, fracture)
%   [G,F,fracture] = gridFracture2D(G, fracture, 'pn1', pv1)
%
% REQUIRED PARAMETERS:
%
%   G  - Grid data structure.
%
%   fl - fracture lines represented by its end points as [x1 y1 x2 y2]. fl
%        will have 1 row per fracture line.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   verbose - Enable output.  Default value dependent upon global verbose
%             settings of function 'mrstVerbose'.
%
% RETURNS:
%
%   G        - Grid data structure with an added cell list G.cells.fracture
%              containing indicator calues for cells containing fractures
%              as well as the line and network indices of those fractures.
%
%   fracture - Structure with the following sub-structures:
%             (a) lines - 1-by-rows(fl) structure with fields network
%                         (network to which each line belongs) and endp
%                         (endpoints of each line as supplied by 'fl')
%             (b) network - stores indices for fracture lines contained in
%                           each network.
%
% NOTE: 
%   This function calls getIndepNetwork and markcells internally.
%
% SEE ALSO:
%   getIndepNetwork, markcells2D

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


opt = struct('verbose', true);
opt = merge_options(opt, varargin{:});

% Get independant fracture lines
dispif(opt.verbose, '\nExtracting independant fracture networks from fracture lines...\n\n');
fracture = getIndepNetwork(fl);

% Mark cells
dispif(opt.verbose, 'Marking cells containing fractures...\n\n');
[G,fracture] = markcells2D(G,fracture); % Designed for 2D as of now


return