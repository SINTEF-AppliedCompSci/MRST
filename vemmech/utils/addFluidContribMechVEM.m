function fbc = addFluidContribMechVEM(G, bc, rock, isdirdofs)
%
%
% SYNOPSIS:
%   function fbc = addFluidContribMechVEM(G, bc, rock, isdirdofs)
%
% DESCRIPTION: Setup the force boundary condition for the mechanics
% corresponding to the fluid boundary condition.
%
% PARAMETERS:
%   G         - Grid structure
%   bc        - Fluid boundary conditions.
%   rock      - Rock structure
%   isdirdofs - Degrees of freedom where in fact Dirichlet boundary
%               conditions are imposed
%
% RETURNS:
%   fbc - Force load for the mechanical system.
%
% EXAMPLE:
%
% SEE ALSO:
%

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


    if (~isempty(bc))
        faces    = bc.face;
        lnn      = G.faces.nodePos(faces + 1) - G.faces.nodePos(faces);
        inodes   = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces + 1) - 1);
        nodes    = G.faces.nodes(inodes);
        lfacenum = rldecode([1:numel(faces)]', lnn);%#ok
        facenum  = rldecode(faces, lnn);
        
        %assert boundary face
        assert(all(sum(G.faces.neighbors(faces, :) == 0, 2) == 1));
        lcells = sum(G.faces.neighbors(faces, :), 2);
        N = bsxfun(@times, G.faces.normals(faces, :), (2*(G.faces.neighbors(faces, 1) == lcells)-1));

        if (G.griddim == 2)
            qf = N/2;
            qf = qf(lfacenum, :);
        else
            % here one should weith the normal with teh subface area
            % this could have been taken from precomputed weights??
            [qc qfs] = calculateQC(G);%#ok
            N        = bsxfun(@rdivide, N, G.faces.areas(faces));
            qf       = bsxfun(@times, N(lfacenum, :), qfs(inodes));
            %error()
        end
        
        % find nodes corresponding to forces
        dofs    = mcolon(G.griddim*(nodes-1)+1, G.griddim*(nodes-1)+G.griddim);
        assert(all(any(G.faces.neighbors(facenum, :) == 0, 2)));
        lfalpha = rock.alpha(sum(G.faces.neighbors(facenum, :), 2));%scale the face with alpha
        force   = bsxfun(@times, qf, lfalpha.*rldecode(bc.value, lnn));
        % transform to dofs numbering
        force = force';
        ndof  = G.nodes.num*G.griddim;
        fbc   = accumarray(dofs', force(:)', [ndof, 1]);
        % maps to current degees of freedom
        fbc   = fbc(~isdirdofs);
    else
        ndof  = G.nodes.num*G.griddim;
        fbc   = zeros(ndof, 1);
        fbc   = fbc(~isdirdofs);
    end
end
