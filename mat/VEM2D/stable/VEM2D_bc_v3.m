function [A, b] = VEM2D_bc_v3(G, A, b, bc, k)

NK = G.nodes.num + G.faces.num*(k-1) + G.cells.num*k*(k-1)/2;

edges   = bc.face(strcmp(bc.type,'flux'));
edgeLengths = G.faces.areas(edges);

nodeNum = mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges +1)-1);
nodes   = G.faces.nodes(nodeNum);
nN = numel(nodes);
uNodes = unique(nodes);
nUN = numel(uNodes);
S = (repmat(nodes,1,nUN) == repmat(uNodes',nN,1))';

vals = bc.value(strcmp(bc.type,'flux'),:);

if k == 1
    vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths/3];
    vals = bsxfun(@plus, vals(:,1:2), vals(:,3));
    vals = reshape(vals',[],1);
    vals = S*vals;
    dofVec = uNodes';
elseif k == 2
    vals = vals.*[edgeLengths/6, edgeLengths/6, edgeLengths*2/3];
    vals = [S*reshape(vals(:,1:2)',[],1); vals(:,3)];
    dofVec = [uNodes', edges' + G.nodes.num];
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
    