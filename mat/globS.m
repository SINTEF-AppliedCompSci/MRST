function S = globS(G, nx, ny)

    Nc = G.cells.num;
    Ne = G.faces.num;
    Nn = G.nodes.num;
    Ndof = Nn + Ne + Nc;
    f = -1;
    
    S = sparse(Ndof,Ndof);
    
    for c = 1:Nc;
        
        nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;
        nodes = G.cells.nodes(nodeNum);
        X = G.nodes.coords(nodes,:);
                                %   Find edge midpoints.
        faceNum = G.cells.facePos(c) : G.cells.facePos(c+1)-1;
        faces = G.cells.faces(faceNum);
        Xmid = G.faces.centroids(faces,:);
                                %   Find boundary edges.
        neighbors = G.faces.neighbors(faces,:);
        boundaryEdges = sum(neighbors == 0, 2);
        %be = boundaryEdges(faces)
                                %   Find edge normals, fix orientation.
        normals = G.faces.normals(faces,:);
        m = (-ones(length(normals),1)).^(neighbors(:,1) ~= c);
        normals = [m,m].*normals;
        edgeLengths = G.faces.areas;
                                %   Find volume.
        vol = G.cells.volumes(c);
        hK = sqrt(2);           %   FIX
        
        [Sl, bl] = locS(X, Xmid, edgeLengths, normals, boundaryEdges, vol, f, hK);

    
        %  MAP LOCAL STIFFNESS MATRIX TO GLOBAL STIFFNESS MATRIX            %%
        
        dofVec = [nodes', faces + Nn, Nn + Ne + c];
        S(dofVec, dofVec) = S(dofVec, dofVec) + Sl;

    end
