function G = VEM2D_makeInternalBoundary(G, faces)
%   Makes internal boundary in gird G, compatible with VEM2D solver.
%
%   SYNOPSIS:
%       G = VEM2D_makeInternalBoundary(G, faces)
%
%   DESCRIPTION:
%       Makes internal boundary in gird G, compatible with VEM2D solver, by
%       duplicating faces and nodes on specified internal boundary faces,
%       and updating maps from faces to nodes, cells to faces and cells to
%       nodes. Function only supports one boundary at a time.
%
%   REQUIRED PARAMETERS:
%       G       - MRST grid with sorted edges, G = sortEdges(G), and VEM2D
%                 geometry, G = computeVEM2DGEometry(G).
%       faces   - faces making up internal boundary.

%   RETURNS:
%       G       - Grid updated with internal boundary.
%-----------------------------------------------------------------ØSK-2016-

%{
   Copyright (C) 2016 Øystein Strengehagen Klemetsdal. See Copyright.txt
   for details.
%}

faces = faces(G.faces.neighbors(faces,1) ~= 0 & ...
              G.faces.neighbors(faces,2) ~=0);
faces = unique(faces);
nE = numel(faces);

neighbors = G.faces.neighbors(faces,:);
G.faces.neighbors(faces,1) = 0;

nodeNum = mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1);
nodes = G.faces.nodes(nodeNum);

nN = numel(nodes);
[uNodes, IU, ID] = unique(nodes);
newCoords = G.nodes.coords(nodes(IU),:);
nUN = numel(uNodes);
%ii = mcolon(2:2:nUN,1:2:nUN,-1);
G.nodes.coords = [G.nodes.coords; newCoords];
newNodes = ((1:nUN) + G.nodes.num)';
G.nodes.num = G.nodes.num + nUN;

cells = neighbors(:,1);
neighbors(:,2) = 0;
G.faces.neighbors = [G.faces.neighbors; neighbors];
G.faces.nodes = [G.faces.nodes; newNodes(ID)];
newNodePos = (G.faces.nodePos(end)+2:2:G.faces.nodePos(end)+nN)';
G.faces.nodePos = [G.faces.nodePos; newNodePos];
% G.faces.tag = [G.faces.tag; zeros(nUN,1)];
newFaces = ((1:nE) + G.faces.num)';
G.faces.num = G.faces.num + nE;


faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
cellFaces = G.cells.faces(faceNum);
cellFaces(ismember(cellFaces, faces)) = newFaces';
G.cells.faces(faceNum) = cellFaces;

nodeNum = mcolon(G.faces.nodePos(cellFaces),G.faces.nodePos(cellFaces+1)-1);
cellFaceNodes = G.faces.nodes(nodeNum);
[cNodes,loc] = ismember(cellFaceNodes,uNodes);

cellFaceNodes(cNodes) = newNodes(loc(cNodes))';
G.faces.nodes(nodeNum) = cellFaceNodes;

% for i = 1:numel(cellFaceNodes)
%     plot(G.nodes.coords(nodes,1),G.nodes.coords(nodes,2),'.','markersize',20)
% end

nodeNum = mcolon(G.cells.nodePos(cells),G.cells.nodePos(cells+1)-1);
cellNodes = G.cells.nodes(nodeNum);
cellNodes(ismember(cellNodes, uNodes)) = newNodes(ID)';
G.cells.nodes(nodeNum) = cellNodes;

end