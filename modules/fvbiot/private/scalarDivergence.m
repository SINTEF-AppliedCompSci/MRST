function div = scalarDivergence(G)
% For a grid G, compute divergence operator (for scalar equations).
%
% G should be an MRST grid, see http://www.sintef.no/projectweb/mrst/ for
% more information.
%
% RETURNS: div - num_cells x num_faces sparse matrix, divergence mapping
%   from face to cell variables, accounting for sign of normal vectors.
%
%{
Parital copyright 2009-2016 SINTEF ICT, Applied Mathematics.
Partial copyright 2016, University of Bergen.

This file is part of FVBiot.

FVBiot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

FVBiot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
cellNo = rldecode(1:G.cells.num, double(diff(G.cells.facePos)), 2).';
cf     = G.cells.faces(:,1);
nc     = G.cells.num;
nf     = G.faces.num;
sgn = 2 * (G.faces.neighbors(cf,1)==cellNo)-1;
div = sparse(cf,cellNo,sgn,nf,nc)';