function bc = pside(bc, G, side, T, varargin)
%Impose temprature boundary condition on global side.
%
% SYNOPSIS:
%   bc = pside(bc, G, side, T)
%   bc = pside(bc, G, side, T, 'pn', pv)
%   bc = pside(bc, G, side, T, I1, I2)
%   bc = pside(bc, G, side, T, I1, I2, 'pn', pv)
%
% PARAMETERS:
%   bc     - Boundary condition structure as defined by function 'addBC'.
%
%   G      - Grid structure as described by grid_structure.  Currently
%            restricted to grids produced by functions cartGrid and
%            tensorGrid.
%
%   side   - Global side from which to extract face indices.  String.  Must
%            (case insensitively) match one of six alias groups:
%
%               1) {'West' , 'XMin', 'Left'  }
%               2) {'East' , 'XMax', 'Right' }
%               3) {'South', 'YMin', 'Back'  }
%               4) {'North', 'YMax', 'Front' }
%               5) {'Upper', 'ZMin', 'Top'   }
%               6) {'Lower', 'ZMax', 'Bottom'}
%
%            These groups correspond to the cardinal directions mentioned
%            as the first alternative in each group.
%
%   T      - temprature value, in units of kelvin, to be applied to the face.
%            Either a scalar or a vector of NUMEL(I1)*NUMEL(I2) values.
%
%   I1, I2 - Cell index ranges for local (in-plane) axes one and two,
%            respectively.  An empty index range ([]) is interpreted as
%            covering the entire corresponding local axis of 'side' in the
%            grid 'G'.  The local axes on a 'side' in 'G' are ordered
%            according to 'X' before 'Y', and 'Y' before 'Z'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   sat    - Fluid composition of fluid injected across inflow faces.
%            An n-by-m array of fluid compositions with 'n' being the
%            number of individual faces specified by (I1,I2) (i.e.,
%            n==NUMEL(I1)*NUMEL(I2)) or one.  If m=3, the columns of 'sat'
%            are interpreted as 1 <-> Aqua, 2 <-> Liquid, 3 <-> Vapor.
%
%            This field is for the benefit of transport solvers such as
%            'blackoilUpwFE' and will be ignored for outflow faces.
%
%            Default value: sat = [] (assume single-phase flow).
%
%   range  - Restricts the search for outer faces to a subset of the cells
%            in the direction perpendicular to that of the face. Example:
%            if side='LEFT', one will only search for outer faces in the
%            cells with logical indexes [range,:,:].
%            Default value: range = [] (do not restrict search).
%
% RETURNS:
%   bc     - Updated boundary condition structure.
%
% EXAMPLE:
%   See simpleBC, simpleSRCandBC.
%
% SEE ALSO:
%   fluxside, addBC, solveIncompFlow, grid_structure.

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

% fake pressure to reuse code
bc = pside(bc, G, side, T, varargin{:});
for i=1:numel(bc.type)
   bc.type{i}='temperature'; 
end
bc=rmfield(bc,'sat');



