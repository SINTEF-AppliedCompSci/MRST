function [intFun, w, x] = makeCellIntegrator(G, cells, degree, type)

    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));

    switch type
            
        case 'tri'
           
            [xr, w, nq] = getQuadratureRule(degree, type);
            b2c = mapBaryToCart(G, cells);
                 
            x = b2c(xr);
            
            vol = getTvolumes(G, cells);
            w = w'.*vol';
            
            ncf = diff(G.cells.facePos);
            
            
            [ii, jj] = blockDiagIndex(ones(numel(cells), 1), ncf(cells)*nq);
            W = sparse(ii, jj, w(:));
            
            intFun = @(fun) W*fun(x);
            
        case 'div'
    
            nodes = reshape(nodes, 2, [])';

            x0 = G.nodes.coords(nodes(:,1),:);
            x1 = G.nodes.coords(nodes(:,2),:);

            [xr, w] = getQuadratureRule(degree);
            xr = reshape(xr, 1, 1, []);
            w = reshape(w, 1, []);

            x = ((x1 - x0).*xr +  (x1 + x0))/2;

            ncf = diff(G.cells.facePos);
            sign = 1 - 2*(G.faces.neighbors(faces,1) ~= rldecode(cells, ncf(cells), 1));
            nx = G.faces.normals(faces,1).*sign/2;

            val = @(fun) integrate(fun, cells, faces, w.*nx, x, ncf);
            
    end
                  
            
end

function val = integrate(fun, cells, faces, w, x, ncf)

    val = zeros(numel(faces), numel(w));
    for wNo = 1:numel(w)
        val(:,wNo) = fun(x(:,:,wNo));
    end
    
    val = sum(val.*w, 2);
    val = accumarray(rldecode((1:numel(cells))',ncf(cells),1), val);
    
end

function vol = getTvolumes(G, cells)

    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
    nodes = reshape(nodes, 2, [])';
    
    ncf = diff(G.cells.facePos);
    xc = rldecode(G.cells.centroids(cells), ncf(cells), 1);
    
%     xc = rldecode(G.cells.centroids, accumarray(ind, nfn(faces)), 1);
    
    v1 = G.nodes.coords(nodes(:,1),:) - xc;
    v2 = G.nodes.coords(nodes(:,2),:) - xc;
    
    vol = abs(v1(:,1).*v2(:,2) - v1(:,2).*v2(:,1))/2;

end