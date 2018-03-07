function G = computeCellDimensions(G)

    n = G.cells.nodes;
    x = G.nodes.coords(n);
    
    ncn = diff(G.cells.nodePos);
    
    cells = rldecode((1:G.cells.num)', ncn, 1);
    nre = G.cells.nodes(mcolon(G.cells.nodePos(cells), G.cells.nodePos(cells+1)-1));
    
    nrl = G.cells.nodes;
    nrl = rldecode(nrl, rldecode(ncn, ncn, 1), 1);
    
    xrl = G.nodes.coords(nrl,:);
    xre = G.nodes.coords(nre,:);
    
    ii = rldecode((1:G.cells.num)', ncn.^2, 1);
    jj = mcolon(1, ncn.^2)';
    
%     [ii, jj] = blockDiagIndex(ncn.^2, ones(G.cells.num, 1));
%     jj = mod(jj-1, max(ncn.^2));
    
    dx = zeros(G.cells.num, G.griddim);
    for dNo = 1:G.griddim
    
        dX = sparse(ii, jj, abs(xrl(:, dNo) - xre(:, dNo)));
        dx(:, dNo) = max(dX, [], 2);
        
    end
    
    G.cells.dx = dx;
   
%     for dNo = 1:G.griddim
%         
%         
%         
%         X = sparse(ii, jj, x(:,dNo));
%         
%         
%     end
%         
    


%     f    = G.cells.faces(:,1);
%     n    = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
%     
%     
%     swap = G.faces.neighbors(f,1) ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1);
%     n    = reshape(n, 2, [])';
%     n(swap,:) = n(swap, [2,1]); n = n(:,1);
%     
% %     ncn = diff(G.cells.facePos);
% %     
% %     x = G.nodes.coords(n,:);
% %     
% %     
% %     
% %     ii = mcolon(G.cells.facePos(1:end-1):G.cells.facePos(end)-1);
%     
%     G.cells.dx = zeros(G.cells.num,G.griddim);
% 
%     for c = 1:G.cells.num
%         f = G.cells.faces(G.cells.facePos(c):G.cells.facePos(c+1)-1);
%         n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
%     
%         swap = G.faces.neighbors(f,1) ~= c;
%         n = reshape(n,2,[])';
%         n(swap, :) = n(swap, [2,1]); n = n(:,1);
%         
%         x = G.nodes.coords(n,:);
%         
%         dx = max(max(pdist2(x(:,1), x(:,1), 'euclidean')));
%         dy = max(max(pdist2(x(:,2), x(:,2), 'euclidean')));
%         
%         G.cells.dx(c,:) = [dx, dy];
%         
%     end
end