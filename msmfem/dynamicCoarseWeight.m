function t = dynamicCoarseWeight(cellno, BI, ct, phi, v, dpdt)
%Compute synthetic multiscale weighting function.
%
% SYNOPSIS:
%   theta = dynamicCoarseWeight(cellno, BI, ct, phi, v, dpdt)
%
% PARAMETERS:
%   cellno - A map giving the global cell number of any given half-face.
%            That is, cellno(i) is the global cell to which half-face 'i'
%            is connected.  Typically,
%
%                nc     = G.cells.num;
%                nf     = double(G.cells.numFaces);
%                cellno = rldecode((1 : nc) .', nf);
%
%            when 'G' is the grid of a reservoir model.
%
%   BI     - Inverse mass matrix, correctly updated for effects of total
%            mobility, such that the expression
%
%                f' * (BI \ f)
%
%            gives the (squared) energy norm of a function 'f' represented
%            on all half faces of the model.
%
%   ct     - Total compressibility.  One scalar value for each global cell
%            in the reservoir model.
%
%   phi    - Porosity.  One scalar value for each global cell in the model.
%
%   v      - Current half-face fluxes for all cells in the model.
%
%   dpdt   - An approximate time derivative of the cell pressure values.
%            One scalar value for each cell in the grid.
%
% RETURNS:
%   theta - A \theta function suitable for passing in option pair
%           ('BasisWeighting',theta) to function generateCoarseSystem.
%
% SEE ALSO:
%   `computeMimeticIP`, `generateCoarseSystem`, `pvt`, `rldecode`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


t = phi .* ct .* dpdt;
t = accumarray(cellno, v .* (BI \ v)) - t;
