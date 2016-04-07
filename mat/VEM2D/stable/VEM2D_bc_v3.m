function [A, b] = VEM2D_bc_v3(G, A, b, bc, k)

NK = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

edges   = bc.face(strcmp(bc.type,'flux'));
nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
nN = size(nodes,1);
edgeLengths = G.faces.areas(edges);
vals = bc.value(strcmp(bc.type,'flux'),:);
nodes = reshape(nodes,2,[])';
[~,I1] = sort(nodes(:,1));
[nodes,I2] = sort(nodes(:,2));
if k == 1
vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths/3];
    vals = bsxfun(@plus, vals(:,1:2), vals(:,3));
    vals = vals(I1,1) + vals(I2,2);
    dofVec = nodes';
elseif k == 2
    vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths*2/3];
    vals = [vals(I1,1) + vals(I2,2); vals(:,3)];
    dofVec = [nodes', edges' + G.nodes.num];
end
b(dofVec) = b(dofVec) + vals;

edges   = bc.face(strcmp(bc.type,'pressure'));
nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
vals = bc.value(strcmp(bc.type,'pressure'),:);

if k == 1
    vals = reshape(vals(:,1:2)',[],1);
    dofVec = nodes';
elseif k == 2
    dofVec = [nodes', edges' + G.nodes.num];
    vals = [reshape(vals(:,1:2)',[],1); vals(:,3)];
end
b(dofVec) = vals;
I = spdiags(ones(NK,1),0,NK,NK);
A(dofVec,:) = I(dofVec,:);

end

% cellEdges = rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
% cellEdges = cellEdges(edges);
% edgeSign = (-ones(size(edges,1),1)).^(G.faces.neighbors(edges,1) ~= cellEdges);
% 
% nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
% nodes   = G.faces.nodes(nodeNum);
% nodes = reshape(nodes,2,[])';
% nN = size(nodes,1);
% nodes(edgeSign == -1,:) = nodes(edgeSign == -1,2:-1:1);
% nodes = nodes(:,1);
    