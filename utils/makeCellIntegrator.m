function [x, w, nq, ii, jj, cellNo, faceNo] = makeCellIntegrator(G, cells, degree, type)

    switch type
        case 'volume'

            [xr, w, nq] = getQuadratureRule(degree, G.griddim);
            
            if degree <= 1

                x = xr + G.cells.centroids(cells,:);
                [ii, jj] = deal((1:numel(cells))');
                cellNo = cells;
                vol = G.cells.volumes(cells);
                nq = repmat(nq, G.cells.num, 1);

            else

                [b2c, nct, vol] = mapBaryToCart_face(G, cells, 'volume');
                x = b2c(xr);
                
                nq = nq*nct;

                [ii, jj] = blockDiagIndex(ones(numel(cells), 1), nq);
                cellNo = rldecode(cells, nq, 1);

            end

            w = reshape(w.*vol', [], 1);
            faceNo = [];
            
        case 'surface'
            
            [xr, w, nq] = getQuadratureRule(degree, G.griddim-1);
            
            faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
            [bf, bc] = boundaryFaces(G);
            ncbf = sum((1:G.cells.num)' == bc',2);
            ncf = diff(G.cells.facePos);
            faces = faces(~ismember(faces, bf));
            
            if G.griddim == 2
            
                l2c = mapLineToCart(G, faces);
                x = l2c(xr);

                w  = repmat(w , numel(faces), 1);
                w = w/2;
                nct = 1;

                
            else
                
%                 ncf = ncf(cells);
                
                if degree <= 1

                    x = G.faces.centroids(faces,:);
%                     [ii, jj] = deal((1:numel(faces))');
%                     cellNo = rldecode(cells, ncf, 1);
%                     faceNo = faces;
%                     vol = G.faces.areas(faces);
%                     nq = repmat(nq, G.cells.num, 1);
                      nct = 1;

                else

%                     [b2c, vol, nct] = mapBaryToCart(G, cells);
                    [b2c, nct, vol] = mapBaryToCart_face(G, cells, 'surface');
                    x = b2c(xr);

%                     nq = nq*nct;

%                     [ii, jj] = blockDiagIndex(ones(numel(cells), 1), nq);
%                     cellNo = rldecode(cells, nq, 1);

                end
                
%                 w = repmat(w , numel(faces), 1);
                w = repmat(w , sum(nct), 1);
%                 w = w/4;
                
            end
            
            swap = G.faces.neighbors(faces,1) ~= rldecode(cells, ncf(cells) - ncbf(cells), 1);

            sign = 1 - 2*swap;

            sign = rldecode(sign, nct*nq, 1);
            vol = reshape(repmat(vol', nq, 1), [], 1);
%             vol  = rldecode(vol, nct*nq, 1);
%             sign = reshape(repmat(sign', nq, 1), [], 1);
            
%             w = repmat(w , numel(faces), 1);
            %     sign = reshape(repmat(sign'.*G.faces.areas(faces)', nq, 1), [], 1);
            w = sign.*vol.*w;
                
            if size(faces,1) == 1, faces = faces'; end
%             faceNo = reshape(repmat(faces', nq, 1), [], 1);
            faceNo = rldecode(faces, nct*nq, 1);

            
            if G.griddim == 3 && degree > 1
                
                vv = rldecode((1:numel(cells))', ncf(cells) - ncbf(cells), 1);
%                 vv(cells)
%                 vv = rldecode((1:G.cells.num)', ncf(cells) - ncbf(cells), 1);
%                 vv(cells)
                nq = accumarray(vv, nct*nq);

            else
                
                nq = (ncf(cells) - ncbf(cells))*nq;
                
            end
            
            
            [ii, jj] = blockDiagIndex(ones(numel(cells), 1), nq);

            cellNo = rldecode(cells, nq, 1);

            
    end

    
            
end

function [x, w, nq, ii, jj, cellNo, faceNo] = makeCellSurfaceIntegrator(G, cells, degree)

%     opt = struct('exclude_boundary', true);
%     opt = merge_options(opt, varargin{:});
    
    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
    if 1%opt.exclude_boundary
        [bf, bc] = boundaryFaces(G);
        ncbf = sum((1:G.cells.num)' == bc',2);
        faces = faces(~ismember(faces, bf));
    else
        ncbf = zeros(G.cells.num, 1);
    end
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1));
    
    ncf = diff(G.cells.facePos);
    swap = G.faces.neighbors(faces,1) ~= rldecode(cells, ncf(cells) - ncbf(cells), 1);
    
    nodes = reshape(nodes, 2, [])';
%     nodes(swap,:) = nodes(swap, [2,1]);

    x0 = G.nodes.coords(nodes(:,1),:);
    x1 = G.nodes.coords(nodes(:,2),:);

    [xr, w, nq] = getQuadratureRule(degree, G.griddim-1);

    x1 = reshape(repmat(x1', nq, 1), G.griddim, [])';
    x0 = reshape(repmat(x0', nq, 1), G.griddim, [])';
    
    xr = repmat(xr, numel(faces), 1);
    w  = repmat(w , numel(faces), 1);
    
    x = ((x1 - x0).*xr +  (x1 + x0))/2;

    
    sign = 1 - 2*swap;

    sign = reshape(repmat(sign', nq, 1), [], 1);
%     sign = reshape(repmat(sign'.*G.faces.areas(faces)', nq, 1), [], 1);
    w = sign.*w/2;
    
%     faceNo = rldecode(faces, 
    if size(faces,1) == 1, faces = faces'; end
    faceNo = reshape(repmat(faces', nq, 1), [], 1);
    
    nq = (ncf(cells) - ncbf(cells))*nq;
    [ii, jj] = blockDiagIndex(ones(numel(cells), 1), nq);

    cellNo = rldecode(cells, nq, 1);

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