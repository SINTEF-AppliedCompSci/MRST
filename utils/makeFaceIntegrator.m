function [x, w, nq, ii, jj, cellNo, faceNo] = makeFaceIntegrator(G, cells, degree, varargin)

    opt = struct('exclude_boundary', true);
    opt = merge_options(opt, varargin{:});
    
    faces = G.cells.faces(mcolon(G.cells.facePos(cells), G.cells.facePos(cells+1)-1));
    if opt.exclude_boundary
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

%     nx = G.faces.normals(faces,1).*sign;
%     nx = reshape(repmat(nx', nq, 1), 1, [])';
%     fa = reshape(repmat((sign.*G.faces.areas(faces))', nq, 1), [], 1);
%     w = fa.*w/2;
%     w  = w.*nx/2;