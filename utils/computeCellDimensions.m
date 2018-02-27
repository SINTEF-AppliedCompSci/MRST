function G = computeCellDimensions(G)

    
    f = G.cells.faces(:,1);
    n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
    swap = G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
    n = reshape(n, 2, [])';
    n(swap,:) = n(swap, [2,1]); n = n(:,1);
    
%     ncn = diff(G.cells.facePos);
%     
%     x = G.nodes.coords(n,:);
%     
%     
%     
%     ii = mcolon(G.cells.facePos(1:end-1):G.cells.facePos(end)-1);
    
    G.cells.dx = zeros(G.cells.num,G.griddim);

    for c = 1:G.cells.num
        f = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
        n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
    
        swap = G.faces.neighbors(f,1) ~= c;
        n = reshape(n,2,[])';
        n(swap, :) = n(swap, [2,1]); n = n(:,1);
        
        x = G.nodes.coords(n,:);
        
        dx = max(max(pdist2(x(:,1), x(:,1), 'euclidean')));
        dy = max(max(pdist2(x(:,2), x(:,2), 'euclidean')));
        
        G.cells.dx(c,:) = [dx, dy];
        
    end
end