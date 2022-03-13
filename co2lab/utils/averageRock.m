function rock2D = averageRock(rock, g_top)
% Average version of rock for use in vertical averaging
%
% SYNOPSIS:
%   rock = averageRock(rock, g_top)
%
% PARAMETERS:
%   rock    - rock structure for 3D grid.  
%
%   g_top   - top surface 2D grid as defined by function 'topSurfaceGrid'.
%
% RETURNS:
%   rock2D  - rock structure with porosities and (lateral) permeability
%             averaged for each column, as well as net-to-gross (ntg)
%             values averaged for each column if ntg is a field of rock
%             structure.

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

c = g_top.columns;
cP = g_top.cells.columnPos;

colNo = rldecode(1:g_top.cells.num, cP(2:end)- cP(1:(end-1)),2).';
pNo = size(rock.perm,2);
if pNo ~= 1
    dispif(mrstVerbose, 'Anisotropic permeability not implemented for VE models, using x-direction.\n')
    rock.perm = rock.perm(:,1);
    pNo = 1;
end
rock2D.poro = accumarray(colNo, rock.poro(c.cells).*c.dz)./accumarray(colNo,c.dz);
rock2D.perm = zeros(g_top.cells.num, pNo);
for i = 1:pNo
    rock2D.perm(:,i) = accumarray(colNo, rock.perm(c.cells,i).*c.dz)./accumarray(colNo,c.dz);
end

%{
X = sparse(colNo, 1:numel(c.cells), 1) * ...
    [rock.poro(c.cells), rock.perm(c.cells,:), ones([numel(c.cells), 1])];

X = bsxfun(@rdivide, X(:,1:end-1), X(:,end));

rock2D.poro = X(:,1);
rock2D.perm = X(:,2:end);
%}


% if net-to-gross is present in rock structure, perform vertical averaging:
if isfield(rock,'ntg')
    rock2D.ntg = accumarray(colNo, rock.ntg(c.cells).*c.dz)./accumarray(colNo,c.dz);
end


end
