function G = computeVEMGeometry(G)
%Computes VEM geometry of MRST grid G.
%
%   SYNOPSIS:
%       G = computeVEMGeometry(G)
%
%   DESCRIPTION:
%       Computes geometry using MRST functions G = computeGeometry(G) and G
%       = mrstGridWithFullMappings(G), and computes edge data and cell
%       diameters. Edge data is organized in the same way as face data in
%       2D.
%
%   REQUIRED PARAMETERS:
%       G   - MRST grid.
%
%   RETURNS:
%       G   - Grid with computed VEM geometry.
%
% SEE ALSO:
%   computeVirtualIP, incompVEM.

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

%   Written by Øystein Strengehagen Klemetsdal, SINTEF/NTNU, 2016.
    
    require vemmech
    
    G = computeGeometry(G);
    G = createAugmentedGrid(G);
    
    if G.griddim == 2
        %   In 2D we only need to calculate cell diameters.
        
        x = G.nodes.coords(G.cells.nodes,:);
        ncn = diff(G.cells.nodePos);
        l = rldecode([0;cumsum(ncn(1:end-1))]+1, diff(G.cells.nodePos),1);
        u = rldecode(cumsum(ncn(1:end))+1, diff(G.cells.nodePos),1)-1;
        ii = mcolon(l,u);
        cellDiameters = sqrt(sum((x(ii,:) - rldecode(x, rldecode(ncn, ncn,1))).^2, 2));
        [ii, jj] = blockDiagIndex(ncn.^2, ones(G.cells.num,1));
        cellDiameters = full(max(sparse(ii,jj,cellDiameters),[],1)');
        G.cells.('diameters') = cellDiameters;
        
    else
        %   In 3D, we also need edge data, and face diameters.

        %   Calculate edge lengths and edge centroids
        nodeNum = mcolon(G.edges.nodePos(1:end-1),G.edges.nodePos(2:end)-1);
        nodes = G.edges.nodes(nodeNum);
        edgeVec   = G.nodes.coords(nodes(2:2:end),:) -  ...
                     G.nodes.coords(nodes(1:2:end-1),:);
        lengths   = sqrt(sum(edgeVec.^2,2));
        centroids = (G.nodes.coords(nodes(2:2:end),:) +  ...
                     G.nodes.coords(nodes(1:2:end-1),:))./2;
        clear nodeNum nodes edgeVec

        %   Calculate edge normals, given for each edge of each face
        faceNormals = G.faces.normals;
        edgeNum = mcolon(G.faces.edgePos(1:end-1),G.faces.edgePos(2:end)-1);
        edges = G.faces.edges(edgeNum);
        signs = G.faces.edgeSign(edgeNum);
        nodeNum = mcolon(G.edges.nodePos(edges),G.edges.nodePos(edges+1)-1);
        nodes = G.edges.nodes(nodeNum);
        edgeVec = G.nodes.coords(nodes(2:2:end),:)-G.nodes.coords(nodes(1:2:end-1),:);
        edgeVec = edgeVec.*repmat(signs,1,3);
        normals = cross(edgeVec, rldecode(faceNormals, diff(G.faces.edgePos), 1));
        normals = bsxfun(@times, normals, 1./sqrt(sum(normals.^2,2)));

        clear faceNormals edgeNum edges signs nodeNum nodes edgeVec

        %   Create mapping G.cells.edges

        ii = rldecode((1:G.edges.num)', diff(G.edges.nodePos),1);
        jj = G.edges.nodes;
        E = sparse(ii, jj, 1);

        ii = G.cells.nodes;
        jj = rldecode((1:G.cells.num)', diff(G.cells.nodePos), 1);
        C = sparse(ii,jj,1);

        [cellEdges,c] = find(E*C == 2);

        edgePos = [0; find(diff([c; c(end)+1]) > 0)] + 1;

        clear nodes edgeNodes edges

        x = G.nodes.coords(G.cells.nodes,:);
        ncn = diff(G.cells.nodePos);
        l = rldecode([0;cumsum(ncn(1:end-1))]+1, diff(G.cells.nodePos),1);
        u = rldecode(cumsum(ncn(1:end))+1, diff(G.cells.nodePos),1)-1;
        ii = mcolon(l,u);
        cellDiameters = sqrt(sum((x(ii,:) - rldecode(x, rldecode(ncn, ncn,1))).^2, 2));
        [ii, jj] = blockDiagIndex(ncn.^2, ones(G.cells.num,1));
        cellDiameters = full(max(sparse(ii,jj,cellDiameters),[],1)');

        x = G.nodes.coords(G.faces.nodes,:);
        nfn = diff(G.faces.nodePos);
        l = rldecode([0;cumsum(nfn(1:end-1))]+1, diff(G.faces.nodePos),1);
        u = rldecode(cumsum(nfn(1:end))+1, diff(G.faces.nodePos),1)-1;
        ii = mcolon(l,u);
        faceDiameters = sqrt(sum((x(ii,:) - rldecode(x, rldecode(nfn, nfn,1))).^2, 2));
        [ii, jj] = blockDiagIndex(nfn.^2, ones(G.faces.num,1));
        faceDiameters = full(max(sparse(ii,jj,faceDiameters),[],1)');

        G.edges.('lengths')     = lengths;
        G.edges.('centroids')   = centroids;
        G.faces.('edgeNormals') = normals;
        G.cells.('edges')       = cellEdges;
        G.cells.('edgePos')     = edgePos;
        G.cells.('diameters')   = cellDiameters;
        G.faces.('diameters')   = faceDiameters;
    end

end
