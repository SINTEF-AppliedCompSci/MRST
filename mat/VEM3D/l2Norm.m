function norm = l2Norm(G,e)


nN = G.nodes.num;
nE = G.edges.num;
nF = G.faces.num;
nK = G.cells.num;

nodeNum = mcolon(G.cells.nodePos(1:end-1),G.cells.nodePos(2:end)-1);
nodes = G.cells.nodes(nodeNum);
if size(nodes,1) == 1
    nodes = nodes';
end

edgeNum= mcolon(G.cells.edgePos(1:end-1),G.cells.edgePos(2:end)-1);
edges= G.cells.edges(edgeNum);
if size(edges,1) == 1
    edges= edges';
end

faceNum = mcolon(G.cells.facePos(1:end-1),G.cells.facePos(2:end)-1);
faces = G.cells.faces(faceNum);
if size(faces,1) == 1
    faces= faces';
end


NK = diff(G.cells.nodePos) + diff(G.cells.edgePos) + diff(G.cells.facePos) + 1;

nodeStep = [0; cumsum(diff(G.cells.nodePos))];
edgeStep = [0; cumsum(diff(G.cells.edgePos))];
faceStep =  [0; cumsum(diff(G.cells.facePos))];
cellStep = (0:nK)';

ii = [mcolon(...
      G.cells.nodePos(1:end-1) + edgeStep(1:end-1) + ...
             faceStep(1:end-1) + cellStep(1:end-1) , ...
      G.cells.nodePos(2:end)-1 + edgeStep(1:end-1) + ...
             faceStep(1:end-1) + cellStep(1:end-1)), ...
      mcolon(...
             nodeStep(2:end)   + G.cells.edgePos(1:end-1) + ...
             faceStep(1:end-1) + cellStep(1:end-1), ...
             nodeStep(2:end)   + G.cells.edgePos(2:end)-1 + ...
             faceStep(1:end-1) + cellStep(1:end-1)), ...
      mcolon(...
             nodeStep(2:end)   + edgeStep(2:end) + ...
      G.cells.facePos(1:end-1) + cellStep(1:end-1), ...
             nodeStep(2:end)   + edgeStep(2:end) + ...
      G.cells.facePos(2:end)-1 + cellStep(1:end-1)), ...
      cumsum(NK)'];

eMax = zeros(sum(NK),1);
jj = [nodes', edges' + nN, faces' + nN + nE, (1:nK) + nN + nE + nF];
eMax(ii) = e(jj);

eMax = cellfun(@max, mat2cell(eMax, NK, 1));
vol = G.cells.volumes;

norm = sqrt(sum(vol.*eMax.^2));

end

