function [isNeumann,isDirichlet] = classifyBoundaryFaces(G, bc)
% Classify faces on the boundary as either Neumann or Dirichlet.
%
%{
Copyright 2015-2016, University of Bergen.

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

Nf = G.faces.num;
isBoundary = any(G.faces.neighbors == 0,2);

isDirichlet = false(Nf,1);
isNeumann = false(Nf,1);

% Neumann condition by default
isNeumann(isBoundary) = true;

% Find Dirichlet boundaries
if ~isempty(bc)
    for iter1 = 1 : numel(bc.type)
        if strcmpi(strtrim(bc.type(iter1)),'pressure') + ...
                strcmpi(strtrim(bc.type(iter1)),'Dirichlet') > 0
            
            isDirichlet(bc.face(iter1)) = 1;
            isNeumann(bc.face(iter1)) = 0;
        end
    end
    isNeumann = isBoundary & isNeumann;
    isDirichlet = isBoundary & isDirichlet;
end