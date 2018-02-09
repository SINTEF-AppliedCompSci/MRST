function [b2c, vol] = mapBaryToCart(G, cells)


    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
    xn    = G.nodes.coords(nodes,:);
    
    ncf = diff(G.cells.facePos);
    nfn = diff(G.faces.nodePos);
    
    ind = rldecode(cells, ncf(cells), 1);
    
    xc = rldecode(G.cells.centroids, accumarray(ind, nfn(faces)), 1);
    
    
    
    xc = rldecode(G.cells.centroids(cells,:), ncf(cells), 1);

    
    x = nan(size(xn,1)*3/2, G.griddim);
    x(G.griddim + 1:G.griddim+1:end,:) = xc;
    x(isnan(x(:,1)), :) = xn;
    x = [x, ones(size(x,1),1)];
    
    
    R = zeros(G.griddim + 1, (G.griddim + 1)*numel(faces));
    R(:,1:(G.griddim+1):end) = reshape(x(:,1),3,[]);
    R(:,2:(G.griddim+1):end) = reshape(x(:,2),3,[]);
    R(:,3:(G.griddim+1):end) = reshape(x(:,3),3,[]);
    
%     T = xn - xc;
%     
%     [ii, jj] = blockDiagIndex(ones(numel(faces),1)*G.griddim);
%     T = sparse(ii, jj, T(:));

    
%     b2c = @(x) x*T;
    
    b2c = @(xb) map(R, xb);
    
end

function x = map(R, xb)

     xx = xb*R;
     x = zeros(3*(size(xx,2)/3), 2);
     for dNo = 1:2
        xtmp = xx(:, dNo:3:end);
        x(:,dNo) = xtmp(:);
     end
     
     
end
    