function [A, b] = VEM2D_bc_v3(G, A, b, bc, k)

NK = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

boundaryEdges = bc.face;
values = bc.value;

edges   = boundaryEdges(strcmp(bc.type,'flux'));
nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
nN = size(nodes,1);
edgeLengths = G.faces.areas(edges);
vals = values(strcmp(bc.type,'flux'));
if k == 1
    int = rldecode(vals.*edgeLengths/2, 2*ones(size(edges,1),1),1);
    [nodes,I] = sort(nodes,1);
    int = int(I);
    int = sum(reshape(int,2,[]),1)';
    dofVec = nodes(1:2:end)';
elseif k == 2
    int = [rldecode(vals.*edgeLengths/6, 2*ones(size(edges,1),1),1); ...
           vals.*edgeLengths*2/3];
    [nodes,I] = sort(nodes,1);
    int(1:nN/2) = sum(reshape(int(I),2,[]),1)';
    dofVec = [nodes', edges' + G.nodes.num];
end
b(dofVec) = b(dofVec) + int;

edges   = boundaryEdges(strcmp(bc.type,'pressure'));

cellEdges = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
cellEdges = cellEdges(edges);
edgeSign = (-ones(size(edges,1),1)).^(G.faces.neighbors(edges,1) ~= cellEdges);

nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
nodes = reshape(nodes,2,[])';
nN = size(nodes,1);
nodes(edgeSign == -1,:) = nodes(edgeSign == -1,2:-1:1);
nodes = nodes(:,1);
vals = values(strcmp(bc.type,'pressure'));
if k == 1
    dofVec = nodes';
elseif k == 2
    dofVec = [nodes', edges' + G.nodes.num];
    vals = repmat(vals,2,1);
end
b(dofVec) = vals;
I = spdiags(ones(NK,1),0,NK,NK);
A(dofVec,:) = I(dofVec,:);

end
    