function bc = interpolateBC(G, bc)

    value = zeros(numel(bc.face), 3);
    
    isNeu = strcmp(bc.type, 'flux');
    value(isNeu,:) = .5*repmat(bc.value(isNeu), 1, 3);
    
    isDir = strcmp(bc.type, 'pressure');
    f = bc.face(isDir);
    fSign = (-ones(numel(f),1)).^(G.faces.neighbors(f,2) ~= 0); 
    
    ve = bc.value(isDir);
    
    nodes = G.faces.nodes(mcolon(G.faces.nodePos(f), ...
                                 G.faces.nodePos(f+1)-1));
    if size(nodes,1) == 1; nodes = nodes'; end
    nodes   = reshape(nodes,2,[])';
    nodes(fSign == -1,:) = nodes(fSign == -1, 2:-1:1);
    
    
    M = bsxfun(@eq, repmat(nodes(:,2)', size(nodes,1), 1), nodes(:,1));
    pn = M*(1:numel(f))';
    nn = M'*(1:numel(f))';
        
    d = sqrt(sum((G.faces.centroids(f( pn(pn~= 0)),:)-G.faces.centroids(f(pn~=0),:)).^2,2));
    a = (ve(pn(pn ~= 0))-ve(pn ~= 0))./d;
    fa = G.faces.areas(f);
    
    vn = zeros(numel(nodes)/2 + nnz(pn==0),1);
    
    isInt = false(numel(vn),1);
    isInt(1) = pn(1);
    j = 2;
    for i = 2:numel(pn)
        isInt(j) = any([pn(i),nn(i-1)]~= 0, 2);
        j = j+1;
        if pn(i) == 0 && nn(i-1) == 0
            j = j+1;
        end
    end
    
    vn(isInt) = a.*fa(pn~=0)/2 + ve(pn~=0);
    
    ii = find(pn==0) - (0:nnz(pn==0)-1)';
    jj = find(pn==0) + (0:nnz(pn==0)-1)';
    vn(jj) = a(ii).*fa(pn==0)/2 + ve(pn==0);
    
    ii = find(nn==0) - (1:nnz(nn==0))';
    jj = find(nn==0) + (1:nnz(nn==0))';
    vn(jj) = -a(ii).*fa(nn==0)/2 + ve(nn==0);
    
    vn = reshape(rldecode(vn, isInt+1, 1),2,[])';
    value(isDir,:) = [vn, ve];
    
    bc.value = value;
      
end