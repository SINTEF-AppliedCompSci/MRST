function [x, w, nq, ii, jj, cellNo] = makeCellIntegrator(G, cells, degree)

    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
           
    [xr, w, nq] = getQuadratureRule(degree, G.griddim);
            
    if degree <= 1
                
        x = xr + G.cells.centroids(cells,:);
        [ii, jj] = deal((1:numel(cells))');
        cellNo = cells;
        vol = G.cells.volumes(cells);
        nq = repmat(nq, G.cells.num, 1);
               
    else

        [b2c, vol, nct] = mapBaryToCart(G, cells);
        x = b2c(xr);
%         vol = getTvolumes(G, cells);
%         ncf = diff(G.cells.facePos);

%         nq = nq*ncf(cells);

        nq = nq*nct;
        
        [ii, jj] = blockDiagIndex(ones(numel(cells), 1), nq);
        cellNo = rldecode(cells, nq, 1);
                
    end
    
    w = reshape(w.*vol', [], 1);
            
end



% function vol = getTvolumes(G, cells)
% 
%     faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
%     nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
%     nodes = reshape(nodes, 2, [])';
%     
%     ncf = diff(G.cells.facePos);
%     xc = rldecode(G.cells.centroids(cells,:), ncf(cells), 1);
%     
% %     xc = rldecode(G.cells.centroids, accumarray(ind, nfn(faces)), 1);
%     
%     v1 = G.nodes.coords(nodes(:,1),:) - xc;
%     v2 = G.nodes.coords(nodes(:,2),:) - xc;
%     
%     vol = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))/2;
% 
% end