function S = globS(G, nx, ny)

    Nc = G.cells.num;
    Ne = G.faces.num;
    Nn = G.nodes.num;
    Ndof = Nn + Ne + Nc;
    f = -1;
    
    S = sparse(Ndof,Ndof);
    
    for c = 1:Nc;
        
        %faces = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1, 1)
      
        nodeNum = G.cells.nodePos(c):G.cells.nodePos(c+1)-1;
        
%         nodes = zeros(2*numel(faces),1);
%         for i = 1:numel(faces)
%             nodes(2*i-1:2*i) = G.faces.nodes(G.faces.nodePos(faces(i)):G.faces.nodePos(faces(i)+1)-1);
%         end
%         nodes = unique(nodes);
        
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
                                %   FIX: NORMALS ARE AREA WHEIGHETD
        normals = G.faces.normals(faces,:);
        m = (-ones(length(normals),1)).^(neighbors(:,1) ~= c);
        normals = [m,m].*normals;
        for i = 1:size(normals,1)
            normals(i,:) = normals(i,:)./norm(normals(i,:));
        end
        edgeLengths = G.faces.areas(faces);
                                %   Find volume.
        vol = G.cells.volumes(c);
        hK = 0;
        n = size(X,1);
        for i = 1:n
            hK = max(norm(repmat(X(i,:),n,1)-X),hK);
        end
        
        [Sl, bl] = locS(X, Xmid, edgeLengths, normals, boundaryEdges, vol, f, hK);

    
        %  MAP LOCAL STIFFNESS MATRIX TO GLOBAL STIFFNESS MATRIX            %%
        
        if size(nodes, 2) == 1
            nodes = nodes';
        end
        if size(faces, 2) == 1
            faces = faces';
        end
        dofVec = [nodes, faces + Nn, Nn + Ne + c];
        S(dofVec, dofVec) = S(dofVec, dofVec) + Sl;

    end

end