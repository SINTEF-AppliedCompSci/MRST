function model = setIsotropicPermeabilityFun(model, varargin)
% this is the 'setter' function for the ModelParameter handling
% an isotropic permeability field.

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

    [location, v] = deal(varargin(1:end-1), varargin{end});
    [nc, nd] = size(model.rock.perm);
    assert(nd == 1); % This function only handles isotropic permeabilities
    
    % recompute transmissibilities
    T_hf = compute_half_trans(model.G, model.operators, v .* model.G.cells.H);
    
    cf = model.G.cells.faces(:,1);
    nf = model.G.faces.num;
    % mapping from from cell-face to face
    M = sparse(cf, (1:numel(cf))', 1, nf, numel(cf));

    % updating model
    model.rock.perm = v;
    model.operators.T_all = 1./(M * (1./T_hf));
    model.operators.T = model.operators.T_all(model.operators.internalConn);
end

function T_hf = compute_half_trans(G, operators, perm)

    % Mappings from cells to its faces
    cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    faces = G.cells.faces(:,1);
    % Vector from cell to face centroid
    C = G.faces.centroids(faces,:) - G.cells.centroids(cells,:);
    % Oriented normals
    sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
    N   = bsxfun(@times, sgn, G.faces.normals(faces, :));
    % Make function
    cn  = sum(C.*N,2)./sum(C.*C,2);
    
    T_hf = cn .* perm(cells);
    
end
