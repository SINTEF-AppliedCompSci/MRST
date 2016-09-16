function state = calculateCellPressure(state, G, S)
% 
%     NP = diff(G.cells.nodePos);
%     
%     ncn = diff(G.cells.nodePos);
%     nce = diff(G.cells.edgePos);
%     ncf = diff(G.cells.facePos);
%     
%     NPf = rldecode(NP, ncf, 1);
%     
%     nfn = diff(G.faces.nodePos);
%     nfe = diff(G.faces.edgePos);
%     
%     v1 = reshape(S.faceCoords(:,1), 3, []);
%     v2 = reshape(S.faceCoords(:,2), 3, []);
%     
%     f = G.cells.faces(:,1);
%     ccf = rldecode(G.cells.centroids, diff(G.cells.facePos), 1);
%     
%     cx = trinomialExpansion(v1(1,f)', v2(1,f)', G.faces.centroids(f,1)-ccf(:,1), 1)';
%     cy = trinomialExpansion(v1(2,f)', v2(2,f)', G.faces.centroids(f,2)-ccf(:,2), 1)';
%     cz = trinomialExpansion(v1(3,f)', v2(3,f)', G.faces.centroids(f,3)-ccf(:,3), 1)';
%     
%     PiNstar = squeezeBlockDiag(S.PiNstar, NP, 4, sum(NP));
% 
%     PiNstar = bsxfun(@times, PiNstar, state.nodePressure(G.cells.nodes)');
%     
%     iiN = mcolon(rldecode([0;cumsum(NP(1:end-1))]+1, diff(G.cells.facePos), 1), ...
%                 rldecode(cumsum(NP), diff(G.cells.facePos), 1));
%     PiNstar = PiNstar(:, iiN);
%     
%     fSign = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
%             ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
%     fn = bsxfun(@times, G.faces.normals(f,:), fSign./G.faces.areas(f));
%     PiNstar(2:4,:) = PiNstar(2:4,:).*rldecode(fn',NPf', 2)/2; 
%     
%     alpha = [0 1 0 0]; alpha = alpha +1;
%     PiNstar = bsxfun(@rdivide, PiNstar, alpha');
% 
%     cx = rldecode(cx, NPf', 2);
%     cy = rldecode(cy, NPf', 2);
%     cz = rldecode(cz, NPf', 2);
%     c1 = [PiNstar(1,:); zeros(2,size(PiNstar,2))];
%     c2 = bsxfun(@times, PiNstar(2,:), cx([3,1,2],:));
%     c3 = bsxfun(@times, PiNstar(3,:), cy([3,1,2],:));
%     c4 = bsxfun(@times, PiNstar(4,:), cz([3,1,2],:));
%     c = c1+c2+c3+c4;
%     
%     c = sparseBlockDiag(c, NPf, 2);
%     
%     
%     e  = G.faces.edges;   
%     en = G.faces.edgeNormals;
%     en = bsxfun(@times, en, G.edges.lengths(e));
% 
%     n   = G.edges.nodes(mcolon(G.edges.nodePos(e),G.edges.nodePos(e+1)-1));
%     n   = reshape(n,2,[])';
%     n(G.faces.edgeSign == -1,:) = n(G.faces.edgeSign == -1,2:-1:1);
%     n   = n(:,1);
%     
%     T  = sparseBlockDiag(S.faceCoords, repmat(3,[G.faces.num,1]), 1);
%     
%     x = sparseBlockDiag(G.nodes.coords(n,:)-rldecode(G.faces.centroids, nfn,1) , nfn, 1);    
%     x = squeezeBlockDiag(x*T, nfn, sum(nfn), 2);
%                     
%     ec = sparseBlockDiag(G.edges.centroids(e,:)-rldecode(G.faces.centroids, nfe, 1), nfe, 1);
%     ec = squeezeBlockDiag(ec*T, nfe, sum(nfe), 2);
%     
%     en = sparseBlockDiag(en, nfe, 1);    
%     en = squeezeBlockDiag(en*T, nfe, sum(nfe), 2);
%     enx = en(:,1).*G.edges.lengths(e);
%     
%     alpha = [0 1 0]; beta = [0 0 1];
%     alpha = alpha + 1;
%     
%     mVals = bsxfun(@power, repmat([x(:,1); ec(:,1)],1,3), alpha)...
%           .*bsxfun(@power, repmat([x(:,2); ec(:,2)],1,3), beta);
%     
%     pos = [1;cumsum(nfn)+1];
%     ii = 1:size(x,1); jj = ii;
%     jj(1:end-1) = jj(2:end);
%     jj(cumsum(pos(2:end)-pos(1:end-1))) = ii([1;cumsum(pos(2:end-1) - pos(1:end-2))+1]);
%       
%     mVals = bsxfun(@times, (mVals(ii,:) + mVals(jj,:))/6 + mVals(size(x,1)+1:end,:)*2/3, enx);
%     I = sparseBlockDiag(ones(1, sum(nfn)), nfn, 2); 
%     mVals = I*mVals;
%     int = sparseBlockDiag(mVals(f,:), ones(numel(f), 1), 1)*c;
%     
%     e  = G.faces.edges(mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1));   
%     n   = G.edges.nodes(mcolon(G.edges.nodePos(e),G.edges.nodePos(e+1)-1));
%     n   = reshape(n,2,[])';
%     n(G.faces.edgeSign == -1,:) = n(G.faces.edgeSign == -1,2:-1:1);
%     n   = n(:,1);
%     int = squeezeBlockDiag(int, NPf, 1, sum(NPf))';
%     I = sparseBlockDiag(ones(1, sum(NP.*ncf)), NP.*ncf, 2); 
%     state.cellPressure = full(I*int);
    
    ncn = diff(G.cells.nodePos);
    c = squeezeBlockDiag(S.PiNstar, ncn, 4, sum(ncn));
    c = c(1,:)';
    ii = rldecode((1:G.cells.num)', ncn, 1);
    jj = 1:numel(G.cells.nodes);
    I = sparse(ii,jj,1);
    state.cellPressure = full(I*(c.*state.nodePressure(G.cells.nodes)));
      
end
    
function coeff = trinomialExpansion(a, b, c, n)
    
    if n == 0
        alpha = 0; beta = 0; gamma = 0;
    elseif n == 1
        alpha = [1,0,0]; beta = [0,1,0]; gamma = [0,0,1];
    elseif n == 2
        alpha = [2,1,1,0,0,0]; beta = [0,1,0,2,1,0]; gamma = [0,0,1,0,1,2];
    else
        alpha = [3 2 2 1 1 1 0 0 0 0];
        beta  = [0 1 0 2 0 1 3 2 1 0];
        gamma = [0 0 1 0 2 1 0 1 2 3];
    end
    
    r = size(a,1);     
    coeff = repmat(factorial(n)./(factorial(alpha).*factorial(beta).*factorial(gamma)), r, 1);
    coeff = coeff.*bsxfun(@power, a,repmat(alpha,r,1))...
         .*bsxfun(@power, b,repmat(beta,r,1))...
         .*bsxfun(@power, c,repmat(gamma,r,1));
    
end