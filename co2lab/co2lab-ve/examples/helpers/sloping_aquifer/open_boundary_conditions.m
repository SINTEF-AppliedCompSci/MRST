function bc = open_boundary_conditions(G, pfun, cross_sectional_problem)
%Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    bc = []; % Start with an empty set of boundary faces
    
    if G.griddim == 3
        vface_ind = (G.faces.normals(:,3) == 0); % identify all vertical faces
        depth = G.faces.centroids(:, 3);
    else
        assert(G.griddim == 2); % presumably a VE grid
        vface_ind = true(size(G.faces.normals, 1), 1); % all faces are vertical
        depth = G.faces.z;
    end
    
    bface_ind = (prod(G.faces.neighbors, 2) == 0); % identify all boundary faces 
    bc_face_ix = find(vface_ind & bface_ind); % identify all lateral boundary faces

    if cross_sectional_problem
        remove = G.faces.normals(bc_face_ix, 2) ~= 0 | ...
                 G.faces.centroids(bc_face_ix, 1) == 0;
        bc_face_ix(remove) = [];
    end
    
    
    
    p_face_pressure = pfun(depth(bc_face_ix));
    
    % p_face_pressure = cell_pressures(bc_cell_ix);  % set lateral boundary face
    %                                                % pressure equal to pressure
    %                                                % of corresponding cell

    bc = addBC(bc, ...
               bc_face_ix, 'pressure', ...      % Add hydrostatic pressure         
               p_face_pressure, 'sat', [1, 0]); % conditions to open boundary faces
end
