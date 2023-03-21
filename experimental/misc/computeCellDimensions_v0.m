function G = computeCellDimensions_v0(G)
%Undocumented Utility Function

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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

    ncn = diff(G.cells.nodePos);
    
    cells = rldecode((1:G.cells.num)', ncn, 1);
    nre = G.cells.nodes(mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells+1)-1));
    
    nrl = G.cells.nodes;
    nrl = rldecode(nrl, rldecode(ncn, ncn, 1), 1);
    
    xrl = G.nodes.coords(nrl,:);
    xre = G.nodes.coords(nre,:);
    
    ii = rldecode((1:G.cells.num)', ncn.^2, 1);
    jj = mcolon(1, ncn.^2)';
    
%     [ii, jj] = blockDiagIndex(ncn.^2, ones(G.cells.num, 1));
%     jj = mod(jj-1, max(ncn.^2));
    
    dx = zeros(G.cells.num, G.griddim);
    for dNo = 1:G.griddim
    
        dX = sparse(ii, jj, abs(xrl(:, dNo) - xre(:, dNo)));
        dx(:, dNo) = max(dX, [], 2);
        
    end
    
    G.cells.dx = dx;
    
    %%
    
    if G.griddim == 3
        
    nfn = diff(G.faces.nodePos);
    
    faces = rldecode((1:G.faces.num)', nfn, 1);
    nre = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
    
    nrl = G.faces.nodes;
    nrl = rldecode(nrl, rldecode(nfn, nfn, 1), 1);
    
    xrl = G.nodes.coords(nrl,:);
    xre = G.nodes.coords(nre,:);
    
    
    faceNo = rldecode(faces, nfn(faces));
    xrl = faceCoords(xrl, faceNo, G);
    xre = faceCoords(xre, faceNo, G);
    
    ii = rldecode((1:G.faces.num)', nfn.^2, 1);
    jj = mcolon(1, nfn.^2)';
    
    dx = zeros(G.faces.num, G.griddim-1);
    for dNo = 1:G.griddim-1
    
        dX = sparse(ii, jj, abs(xrl(:, dNo) - xre(:, dNo)));
        dx(:, dNo) = max(dX, [], 2);
        
    end
    
    G.faces.dx = dx;
    
    G = faceCoordSys(G);
    
    end
    
end

function G = faceCoordSys(G)
    
    nfe   = diff(G.faces.edgePos);
    edges = G.faces.edges(cumsum(nfe));
    nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));

    vec1 = G.nodes.coords(nodes(2:2:end),:) - G.nodes.coords(nodes(1:2:end-1),:);
    vec1 = vec1./sqrt(sum(vec1.^2,2));

    n    = G.faces.normals./G.faces.areas;
    vec2 = cross(vec1, n, 2);
    
    G.faces.coordSys = {vec1, vec2};

end


function [x, G] = faceCoords(x, faceNo, G)
            
    nfe   = diff(G.faces.edgePos);
    edges = G.faces.edges(cumsum(nfe));
    edges = edges(faceNo);
    nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));

    vec1 = G.nodes.coords(nodes(2:2:end),:) - G.nodes.coords(nodes(1:2:end-1),:);
    vec1 = vec1./sqrt(sum(vec1.^2,2));

    n    = G.faces.normals(faceNo,:)./G.faces.areas(faceNo);
    vec2 = cross(vec1, n, 2);

    x = [sum(x.*vec1,2), sum(x.*vec2, 2)];
    
    
    
end
