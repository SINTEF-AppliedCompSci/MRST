function bcAD = bc2ADbc(G,bc)
% INTERNAL DEPRECATED FUNCTION: Intentionally undocumented.

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

    bcAD=bc;
    assert(all(any(G.faces.neighbors(bc.face,:)==0,2))); %bc schould be on a boundary
    bc_cell=sum(G.faces.neighbors(bc.face,:),2);    
    if(false)
        % this used full matrix dimensions
        bcAD.cell2bcface=sparse(bc.face,bc_cell,1,G.faces.num,G.cells.num);
        bcAD.bcface2cell=bcAD.cell2bcface';
    else
       % this uses matrixed reduced to size of number of faces, but keep
       % the dimension fore cells. one could use smaller dimesion for cells
       %  but then the matrixes would not be transpose but one had to
       %  introduce a matrix for cells to bc_cells.??
       nbc=numel(bc.face);
       bc_face_index=[1:nbc]';
       bcAD.cell2bcface=sparse(bc_face_index,bc_cell,1,nbc,G.cells.num);
       bcAD.bcface2cell=bcAD.cell2bcface';         
    end
end
