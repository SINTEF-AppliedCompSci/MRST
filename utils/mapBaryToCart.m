function [b2c, nct, vol] = mapBaryToCart(G, cells)

    
    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
    if G.griddim == 2
        edges = faces;
        nodes = G.faces.nodes(mcolon(G.faces.nodePos(edges), G.faces.nodePos(edges+1)-1));
        nct = diff(G.cells.facePos);
        nct = nct(cells);
    else
        edges = G.faces.edges(mcolon(G.faces.edgePos(faces), G.faces.edgePos(faces+1)-1));
        nodes = G.edges.nodes(mcolon(G.edges.nodePos(edges), G.edges.nodePos(edges+1)-1));
        nfe = diff(G.faces.edgePos);
        nct = accumarray(rldecode((1:G.cells.num)', diff(G.cells.facePos), 1), nfe(faces));
        nct = nct(cells);
        
    end
    
    xc = rldecode(G.cells.centroids(cells,:), nct, 1);
    xn = G.nodes.coords(nodes,:);
    
    ntpts = G.griddim+1;
    vec = 1:ntpts:sum(nct)*ntpts-ntpts+1;
    ixn = mcolon(vec, vec+1);
    
    ixc = vec + ntpts - 1;
    
    if G.griddim == 3
        
        xf  = rldecode(G.faces.centroids(faces,:), nfe(faces), 1);
        ixf = vec + ntpts - 2;
        
    else
        
        xf  = [];
        ixf = [];
        
    end
    
    x([ixn, ixf, ixc],:) = [xn; xf; xc];
    
    R = zeros(ntpts, G.griddim*sum(nct));
    for dNo = 1:G.griddim
        R(:,dNo:(G.griddim):end) = reshape(x(:,dNo), ntpts, []);
    end

    b2c = @(xb) map(G, R, xb);
    
    %%
    
    x0 = x(vec,:);
    v = cell(G.griddim,1);
    for dNo = 1:G.griddim
%         ix = dNo:G.griddim:G.griddim*sum(nct(cells)) - G.griddim + dNo;
        v{dNo} = x0 - x(vec + dNo, :);
%         v(ix, :) = 
    end
    
    if G.griddim == 2
        vol = abs(v{1}(:,1).*v{2}(:,2) - v{1}(:,2).*v{2}(:,1))/2;
    else
        vol = abs(dot(cross(v{1}, v{2}, 2), v{3}, 2))/6;
    end
    
end

function x = map(G, R, xb)

     xx = xb*R;
     
     x = zeros(numel(xx)/G.griddim, G.griddim);
     for dNo = 1:G.griddim
        xtmp = xx(:, dNo:G.griddim:end);
        x(:,dNo) = xtmp(:);
     end
          
end