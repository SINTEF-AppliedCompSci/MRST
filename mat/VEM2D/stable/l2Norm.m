function n = l2Norm(G,f,k)

nK = G.cells.num;
n = 0;
for K = 1:nK
    nodeNum = G.cells.nodePos(K):G.cells.nodePos(K+1)-1;
    nodes = G.cells.nodes(nodeNum);
    if size(nodes,1) == 1, nodes = nodes'; end
    if k == 1
        dofVec = nodes;
    else
        faceNum = G.cells.facePos(K):G.cells.facePos(K+1)-1;
        faces = G.cells.faces(faceNum);
        if size(faces,1) == 1, faces = faces'; end
        dofVec = [nodes; faces + G.nodes.num];
    end
    fK = max(abs(f(dofVec)));
    n = n + fK^2*G.cells.volumes(K);
end
n = sqrt(n);
    