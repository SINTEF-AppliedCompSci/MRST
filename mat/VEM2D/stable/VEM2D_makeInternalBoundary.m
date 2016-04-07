function G = VEM2D_makeInternalBoundary(G, faces)


faces = faces(G.faces.neighbors(faces,1) ~= 0 & G.faces.neighbors(faces,2) ~=0);
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
G.nodes.coords = [G.nodes.coords; newCoords];
newNodes = ((1:nUN) + G.nodes.num)';
G.nodes.num = G.nodes.num + nUN;

cells = neighbors(:,1);
neighbors(:,2) = 0;
G.faces.neighbors = [G.faces.neighbors; neighbors];
G.faces.nodes = [G.faces.nodes; newNodes(ID)];
newNodePos = (G.faces.nodePos(end)+2:2:G.faces.nodePos(end)+nN)';
G.faces.nodePos = [G.faces.nodePos; newNodePos];
G.faces.tag = [G.faces.tag; zeros(nUN,1)];
newFaces = ((1:nE) + G.faces.num)';
G.faces.num = G.faces.num + nE;

faceNum = mcolon(G.cells.facePos(cells),G.cells.facePos(cells+1)-1);
cellFaces = G.cells.faces(faceNum);
cellFaces(ismember(cellFaces, faces)) = newFaces';
G.cells.faces(faceNum) = cellFaces;
nodeNum = mcolon(G.cells.nodePos(cells),G.cells.nodePos(cells+1)-1);
cellNodes = G.cells.nodes(nodeNum);
cellNodes(ismember(cellNodes, uNodes)) = newNodes(ID)';
G.cells.nodes(nodeNum) = cellNodes;

G = computeGeometry(G);

% for i = 1:numel(cells)
%     edgeNum = G.cells.facePos(cells(i)):G.cells.facePos(cells(i)+1)-1;
%     edges = G.cells.faces(edgeNum);
%     for j = 1:numel(edges)
%         nodeNum = G.faces.nodePos(edges(j)):G.faces.nodePos(edges(j)+1)-1;
%         nodes = G.faces.nodes(nodeNum);
%         X = G.nodes.coords(nodes,:);
%         plot(X(:,1), X(:,2))
%         hold on
%         pause
%     end
%     hold off
% end

end