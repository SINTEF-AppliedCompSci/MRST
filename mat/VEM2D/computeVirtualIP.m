function S = computeVirtualIP(G,rock,k, varargin)

%%  MERGE INPUT PARAMETRES                                               %%

opt = struct('sigma', []);
opt = merge_options(opt, varargin{:});


%%  CHECK CORRECTNESS OF INPUT

if ~isempty(opt.sigma)
    assert(numel(opt.sigma) == sum(nker));
end

%

%%  CALCULATE 2D PROJECTOION OPERATORS

K = permTensor(rock,G.griddim);
nk   = (k+1)*(k+2)/2;

if G.griddim == 2
    
    %   Calculate projection operators for each cell.
    
    %   Number of nodes and faces for each cell.
    ncn = diff(G.cells.nodePos);
    ncf = diff(G.cells.facePos);
    
    %   Faces for each cell.
    f = G.cells.faces(:,1);
    fn = G.faces.normals(f,:);
    faceSign = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
                  ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
    fn = bsxfun(@times, fn, faceSign);
    if size(f,1) == 1; f = f'; end

    %   Nodes for each face of each cell.
    n = G.faces.nodes(mcolon(G.faces.nodePos(f), ...
                                 G.faces.nodePos(f+1)-1));
    if size(n,1) == 1; n = n'; end
    n   = reshape(n,2,[])';
    n(faceSign == -1,:) = n(faceSign == -1,2:-1:1);
    n   = n(:,1);
    
    %   Function space dimensions.
    
    NP   = ncn + ncf*(k-1) + k*(k-1)/2;
    nker = NP - nk;

    %   Coordinates for degrees of freedom.
    if k == 1
        x = G.nodes.coords(n,:);
    else
        x = [G.nodes.coords(n,:); G.faces.centroids(f,:)];
    end
    
    Kmat = reshape(K', 2, [])';
    
    %   Calculate B and D matrices.
    [B, D, B1, D1] = computeBD2D(G.cells.centroids, G.cells.diameters, ...
                                 G.cells.volumes, ncn, ncf, ...
                                 fn, G.faces.areas(f), G.cells.facePos, ...
                                 x, numel(n), G.cells.nodePos, ...
                                 Kmat, ...
                                 NP, k, G.griddim);
    
    %   Calculate projection operators in monomial (star) and VEM bases.
    M = B*D;
    [ii, jj] = blockDiagIndex(repmat(nk, [G.cells.num ,1]));
    kk = sub2ind(size(M), ii, jj);
    PiNstar = sparse(ii, jj, invv(full(M(kk)), repmat(nk, [G.cells.num, 1])))*B;
    PiN = D*PiNstar;

    clear B D;

    %   Calculate projection operators 
    if k == 1
        PiN1 = (B1*D1)\B1;
    end

    %   Calculate Q matrices-
    Q = zeros(sum(nker.*NP),1);
    PiNPos = [1; cumsum(NP.^2) + 1];
    QPos   = [1; cumsum(NP.*nker)+1];
    ii = blockDiagIndex(NP);
    PiNvec = full(PiN(ii));

    for P = 1:G.cells.num 
        QP = null(reshape(PiNvec(PiNPos(P):PiNPos(P+1)-1), NP(P), NP(P)));
        Q(QPos(P):QPos(P+1)-1) = QP(:);
    end

    [ii,jj] = blockDiagIndex(NP, nker);
    Q = sparse(ii, jj, Q, sum(NP), sum(nker));

    %   Calculate local discrete bilinear forms.
    if isempty(opt.sigma)
        opt.sigma = rldecode(K(:,1)+K(:,4),nker,1);
    end
    sigma = spdiags(opt.sigma, 0, sum(nker), sum(nker));

    M(1:nk:end,:) = 0;
    I = speye(size(PiN,1));
    A = PiNstar'*M*PiNstar + (I-PiN)'*Q*sigma*Q'*(I-PiN);

    %   Make solution struct.
    
    vec = [1, cumsum(NP(1:end-1))' + 1];
    iiN = mcolon(vec, vec + ncn'-1);
    iiF = mcolon(vec + ncn', vec + ncn' + ncf'*(k-1)-1);
    iiP = mcolon(vec + ncn' + ncf', vec + ncn' + ncf' + k*(k-1)/2 -1);
    if k == 1
        dofVec([iiN, iiF, iiP]) = n';
    else
        dofVec([iiN, iiF, iiP]) = [n', f' + G.nodes.num, (1:G.cells.num) + G.nodes.num + G.faces.num*(k-1)];
    end

    S.A = A;
    S.dofVec = dofVec;
    S.PNstar = PiNstar;
    if k == 1
        S.PiN1 = PiN1;
    end
    S.order  = k;
    
else
    
    %%  3D CASE
    
    %%  CALCULATE PROJECTION OPERATORS FOR EACH FACE
    
    nfn = diff(G.faces.nodePos);
    nfe = diff(G.faces.edgePos);
    
    %   Compute local coordinates for each face.

    e  = G.faces.edges;   
    en = G.faces.edgeNormals;
    en = bsxfun(@times, en, G.edges.lengths(e));

    n   = G.edges.nodes(mcolon(G.edges.nodePos(e),G.edges.nodePos(e+1)-1));
    n   = reshape(n,2,[])';
    n(G.faces.edgeSign == -1,:) = n(G.faces.edgeSign == -1,2:-1:1);
    n   = n(:,1);
    nn = numel(n);
    
    x = G.nodes.coords(n,:);

    v1 = (x(G.faces.nodePos(1:end-1)+1,:) - x(G.faces.nodePos(1:end-1),:));
    v1 = bsxfun(@rdivide, v1, sqrt(sum(v1.^2,2)));
    v2 = cross(G.faces.normals,v1,2);
    v2 = bsxfun(@rdivide, v2, sqrt(sum(v2.^2,2)));
    v1 = v1'; v2 = v2';
    T  = sparseBlockDiag([v1(:), v2(:)], repmat(3,[G.faces.num,1]), 1);
    x = sparseBlockDiag(x-rldecode(G.faces.centroids, nfn,1) , nfn, 1);    
    x = squeezeBlockDiag(x*T, nfn, sum(nfn), 2);
                    
    ec = sparseBlockDiag(G.edges.centroids(e,:)-rldecode(G.faces.centroids, nfe, 1), nfe, 1);
    ec = squeezeBlockDiag(ec*T, nfe, sum(nfe), 2);
    
    en = sparseBlockDiag(en, nfe, 1);    
    en = squeezeBlockDiag(en*T, nfe, sum(nfe), 2);
 
    %   Function space dimensions.
    
    
    f = G.cells.faces(:,1);
        
    Kmat = rldecode(K, diff(G.cells.facePos), 1);
    Kmat = Kmat';
    [ii,jj] = blockDiagIndex(repmat(3,[size(Kmat,2), 1]));
    Kmat = sparse(ii, jj, Kmat(:));
    
    T = [v1', v2'];
    T = T(f,:);
    [ii, jj] = blockDiagIndex(repmat(3,[size(T,1),1]),repmat(2,[size(T,1),1]));
    T = T';
    T = sparse(ii, jj, T(:));
    Kmat = squeezeBlockDiag(T'*Kmat*T, repmat(2, [numel(f), 1]), 2*numel(f), 2);
    
%     clear T;
    
    e = G.faces.edges(mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1));
    n = G.faces.nodes(mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1));
   
    
    %   Coordinates for degrees of freedom.
    
    iin = mcolon(G.faces.nodePos(f), G.faces.nodePos(f+1)-1);
    iie = mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1);
    if k == 1
        xx = x(iin,:);
    else
        xx = [x(iin,:); ec(iie,:)];
    end
    
    NP   = nfn + nfe*(k-1) + k*(k-1)/2;
    NP = NP(f);
    
    ePos = diff(G.faces.edgePos);
    ePos = ePos(f);
    ePos = [1;cumsum(ePos)+1];
    nPos = diff(G.faces.nodePos);
    nPos = nPos(f);
    nPos = [1;cumsum(nPos)+1];
    
    
    %   Calculate B and D matrices.
    [BF, DF, ~, ~] = computeBD2D(zeros(numel(f),2), G.faces.diameters(f), ...
                                 G.faces.areas(f), nfn(f), nfe(f), ...
                                 en(iie,:), G.edges.lengths(e), ePos, ...
                                 xx, numel(n), nPos, ...
                                 Kmat, ...
                                 NP, k, G.griddim);

    %   Calculate projection operators in monomial (star) and VEM bases.
    MF = BF*DF;
    [ii, jj] = blockDiagIndex(repmat(nk, [numel(f) ,1]));
    kk = sub2ind(size(MF), ii, jj);
%     PiNFstar = sparse(ii, jj, invv(full(MF(kk)), repmat(nk, [G.cells.num, 1])))*BF;

    PiNFstar = MF\BF;
    
    clear BF DF;
    
    %%  CALCULATE D MATRICES
    
    [m, grad_m, int_m] = retrieveMonomials(3, k);
    ncn = diff(G.cells.nodePos);
    nce = diff(G.cells.edgePos);
    ncf = diff(G.cells.facePos);
    
    NP = diff(G.cells.nodePos) + diff(G.cells.edgePos)*polyDim(k-2,1) ...
       + diff(G.cells.facePos)*polyDim(k-2, 2) + k*(k^2-1)/6*polyDim(k-2,3);
    nk = polyDim(2, G.griddim);
    
    f = G.cells.faces(:,1);
    eNum = mcolon(G.faces.edgePos(f), G.faces.edgePos(f+1)-1);
    e = G.faces.edges(eNum);
    en = G.faces.edgeNormals(eNum,:);
    n = G.edges.nodes(mcolon(G.edges.nodePos(e), G.edges.nodePos(e+1)-1));
    n = reshape(n,2,[])';
    n(G.faces.edgeSign(eNum) == -1,:) = n(G.faces.edgeSign(eNum) == -1,2:-1:1);
    n = n(:,1);
    nn= numel(n);
%         
    if k == 1
        xMon = bsxfun(@rdivide, G.nodes.coords(G.cells.nodes,:) ...
                               - rldecode(G.cells.centroids, ncn,1), ...
                                 rldecode(G.cells.diameters, ncn, 1));
        D = sparseBlockDiag(m(xMon), diff(G.cells.nodePos), 1);
    else
        
        ccf = rldecode(G.cells.centroids, diff(G.cells.facePos), 1);

        cx = trinomialExpansion(v1(1,f)',v2(1,f)', G.faces.centroids(f,1) - ccf(:,1), 1);
        cy = trinomialExpansion(v1(2,f)',v2(2,f)', G.faces.centroids(f,2) - ccf(:,2), 1);
        cz = trinomialExpansion(v1(3,f)',v2(3,f)', G.faces.centroids(f,3) - ccf(:,3), 1);

%         coeff = [coeff, zeros(r, 6)];
%         coeff = [zeros(r, 3), coeff];
        
        alpha = [1 0 0]; beta  = [0 1 0];
        
        [alphaBi, betaBi, c6] = polyProducts(cx, cy, alpha, beta);
        [~      , ~     , c7] = polyProducts(cx, cz, alpha, beta);
        [~      , ~     , c9] = polyProducts(cy, cz, alpha, beta);
        
        alphaBi = alphaBi+1;
        c6 = bsxfun(@rdivide, c6, alphaBi);
        c7 = bsxfun(@rdivide, c7, alphaBi);
        c9 = bsxfun(@rdivide, c9, alphaBi);
        
        c5  = trinomialExpansion(v1(1,f)', v2(1,f)', G.faces.centroids(f,1) - ccf(:,1), 2);
        c8  = trinomialExpansion(v1(2,f)', v2(2,f)', G.faces.centroids(f,2) - ccf(:,2), 2);
        c10 = trinomialExpansion(v1(3,f)', v2(3,f)', G.faces.centroids(f,3) - ccf(:,3), 2);
        
        alphaQuad = [2 1 1 0 0 0];
        betaQuad  = [0 1 0 2 1 0];
        alphaQuad = alphaQuad + 1;
        c5 = bsxfun(@rdivide, c5, alphaQuad);
        c8 = bsxfun(@rdivide, c8, alphaQuad);
        c10 = bsxfun(@rdivide, c10, alphaQuad);
        
        fc = rldecode(G.faces.centroids(f,:), nfn(f), 1);
        x = sparseBlockDiag(G.nodes.coords(n,:)-fc, nfn(f), 1);
        x = squeezeBlockDiag(x*T, nfn(f), sum(nfn(f)), 2);
        ec = sparseBlockDiag(G.edges.centroids(e,:) - fc, nfn(f), 1);
        ec = squeezeBlockDiag(ec*T, nfn(f), sum(nfn(f)), 2);
        en = sparseBlockDiag(en, nfn(f), 1);
        en = squeezeBlockDiag(en*T, nfn(f), sum(nfn(f)), 2);
        enx = en(:,1).*G.edges.lengths(e);
        
        pos = [1;cumsum(nfn(f))+1];
        ii = 1:size(x,1); jj = ii;
        jj(1:end-1) = jj(2:end);
        jj(cumsum(pos(2:end)-pos(1:end-1))) = ii([1;cumsum(pos(2:end-1) - pos(1:end-2))+1]);
        
        mVals = bsxfun(@power, repmat([x(:,1); ec(:,1)],1,6), alphaQuad)...
              .*bsxfun(@power, repmat([x(:,2); ec(:,2)],1,6), betaQuad);
        mVals = bsxfun(@times, (mVals(ii,:) + mVals(jj,:))/6 + mVals(size(x,1)+1:end,:)*2/3, enx);
        mVals = sparseBlockDiag(mVals, nfn(f), 1);
        
        I = sparseBlockDiag(ones(1, sum(nfn(f))), nfn(f), 2); 
        cd = rldecode(G.cells.diameters, diff(G.cells.facePos), 1);
        
        c5 = c5';
        m5fInt = I*mVals*c5(:)./cd.^2;
        
        c8 = c8';
        m8fInt = I*mVals*c8(:)./cd.^2;
        
        c10 = c10';
        m10fInt = I*mVals*c10(:)./cd.^2;
        
        mVals = bsxfun(@power, repmat([x(:,1); ec(:,1)],1,3*3), alphaBi)...
              .*bsxfun(@power, repmat([x(:,2); ec(:,2)],1,3*3), betaBi);
        mVals = bsxfun(@times, (mVals(ii,:) + mVals(jj,:))/6 + mVals(size(x,1)+1:end,:)*2/3, enx);
        mVals = sparseBlockDiag(mVals, nfn(f), 1);
        
        I = sparseBlockDiag(ones(1, sum(nfn(f))), nfn(f), 2); 
        cd = rldecode(G.cells.diameters, diff(G.cells.facePos), 1);
        
        c6 = c6';
        m6fInt = I*mVals*c6(:)./cd.^2;
        
        c7 = c7';
        m7fInt = I*mVals*c7(:)./cd.^2;
        
        c9 = c9';
        m9fInt = I*mVals*c9(:)./cd.^2;
        
        m5fI = []; m8fI = []; m10fI = []; m6fI = []; m7fI = []; m9fI = [];
        
%         for P = 1:G.cells.num
%             
%             ff = G.cells.faces(G.cells.facePos(P):G.cells.facePos(P+1)-1);
%             
%             mm = @(X) ((X(:,1)-G.cells.centroids(P,1))/G.cells.diameters(P)).^2;
%             m5fI = [m5fI; polygonInt3D(G, ff, mm, 2)];
%             
%             mm = @(X) ((X(:,2)-G.cells.centroids(P,2))/G.cells.diameters(P)).^2;
%             m8fI = [m8fI; polygonInt3D(G, ff, mm, 2)];
%             
%             mm = @(X) ((X(:,3)-G.cells.centroids(P,3))/G.cells.diameters(P)).^2;
%             m10fI = [m10fI; polygonInt3D(G, ff, mm, 2)];  
%             
%             mm = @(X) (X(:,1)-G.cells.centroids(P,1)).*(X(:,2)-G.cells.centroids(P,2))...
%                       ./G.cells.diameters(P)^2;
%             m6fI = [m6fI; polygonInt3D(G, ff, mm, 2)];
%             
%             mm = @(X) (X(:,1)-G.cells.centroids(P,1)).*(X(:,3)-G.cells.centroids(P,3))...
%                       ./G.cells.diameters(P)^2;
%             m7fI = [m7fI; polygonInt3D(G, ff, mm, 2)];
%             
%             mm = @(X) (X(:,2)-G.cells.centroids(P,2)).*(X(:,3)-G.cells.centroids(P,3))...
%                       ./G.cells.diameters(P)^2;
%             m9fI = [m9fI; polygonInt3D(G, ff, mm, 2)];
%             
%         end 
%         
%         norm(m5fI - m5fInt)
%         norm(m6fI - m6fInt)
%         norm(m7fI - m7fInt)
%         norm(m8fI - m8fInt)
%         norm(m9fI - m9fInt)
%         norm(m10fI - m10fInt)
        
        
        fSign = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
                  ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
        fn = bsxfun(@times, G.faces.normals(f,:), fSign./G.faces.areas(f));
        
        cx = [cx, zeros(numel(f), polyDim(2, 3) - polyDim(1,3))];
        cy = [cy, zeros(numel(f), polyDim(2, 3) - polyDim(1,3))];
        cz = [cz, zeros(numel(f), polyDim(2, 3) - polyDim(1,3))];
        c5 = [zeros(numel(f), polyDim(1, 3)-1), bsxfun(@times, c5', alphaQuad)];
        c8 = [zeros(numel(f), polyDim(1, 3)-1), bsxfun(@times, c8', alphaQuad)];
        
        alpha = [1 0 0 2 1 1 0 0 0]; beta = [0 1 0 0 1 0 2 1 0];
        
        [alphaBi, betaBi, c6] = polyProducts(c5, cy, alpha, beta);
        [~      , ~     , c7] = polyProducts(c5, cz, alpha, beta);
        [~      , ~     , c9] = polyProducts(c8, cz, alpha, beta);
        
        c6 = bsxfun(@times, c6, fn(:,1)/2);
        c7 = bsxfun(@times, c7, fn(:,1)/2);
        c9 = bsxfun(@times, c9, fn(:,2)/2);
         
        alphaBi = alphaBi +1;
        
        c6 = bsxfun(@rdivide, c6, alphaBi);
        c7 = bsxfun(@rdivide, c7, alphaBi);
        c9 = bsxfun(@rdivide, c9, alphaBi);
        
        c5 = trinomialExpansion(v1(1,f)', v2(1,f)', G.faces.centroids(f,1) - ccf(:,1), 3);
        c8 = trinomialExpansion(v1(2,f)', v2(2,f)', G.faces.centroids(f,2) - ccf(:,2), 3);
        c10 = trinomialExpansion(v1(3,f)', v2(3,f)', G.faces.centroids(f,3) - ccf(:,3), 3);
        
        c5 = bsxfun(@times, c5, fn(:,1)/3);
        c8 = bsxfun(@times, c8, fn(:,2)/3);
        c10 = bsxfun(@times, c10, fn(:,3)/3);
        
        alphaQuad = [3 2 2 1 1 1 0 0 0 0];
        betaQuad  = [0 1 0 2 0 1 3 2 1 0];
        alphaQuad = alphaQuad + 1;
        
        
        c5 = bsxfun(@rdivide, c5, alphaQuad);
        c8 = bsxfun(@rdivide, c8, alphaQuad);
        c10 = bsxfun(@rdivide, c10, alphaQuad);
        
        pos = [1;cumsum(nfn(f))+1];
        ii = 1:size(x,1); jj = ii;
        jj(1:end-1) = jj(2:end);
        jj(cumsum(pos(2:end)-pos(1:end-1))) = ii([1;cumsum(pos(2:end-1) - pos(1:end-2))+1]);
        
        
        eVec = G.nodes.coords(G.edges.nodes(2:2:end),:)...
             - G.nodes.coords(G.edges.nodes(1:2:end),:);
        xq1 = G.edges.centroids - .5*sqrt(1/5)*eVec;
        xq2 = G.edges.centroids + .5*sqrt(1/5)*eVec;
        
        xq1 = sparseBlockDiag(xq1(e,:)-fc, nfn(f), 1);
        xq1 = squeezeBlockDiag(xq1*T, nfn(f), sum(nfn(f)), 2);
        
        xq2 = sparseBlockDiag(xq2(e,:)-fc, nfn(f), 1);
        xq2 = squeezeBlockDiag(xq2*T, nfn(f), sum(nfn(f)), 2);
        
        
        nn = size(x,1);
        mVals = bsxfun(@power, repmat([x(:,1); xq1(:,1); xq2(:,1)],1,10), alphaQuad)...
              .*bsxfun(@power, repmat([x(:,2); xq1(:,2); xq2(:,2)],1,10), betaQuad);
        mVals = bsxfun(@times, (mVals(ii,:) + mVals(jj,:))/12 + (mVals(nn + 1:2*nn,:) + mVals(2*nn+1:3*nn,:))*5/12, enx);
        mVals = sparseBlockDiag(mVals, nfn(f), 1);
        
        If = sparseBlockDiag(ones(1, sum(nfn(f))), nfn(f), 2); 
        Ic = sparseBlockDiag(ones(1, sum(ncf)), ncf, 2); 
        
        c5 = c5';
        m5cInt = Ic*If*mVals*c5(:)./G.cells.diameters.^2;
        
        c8 = c8';
        m8cInt = Ic*If*mVals*c8(:)./G.cells.diameters.^2;
        
        c10 = c10';
        m10cInt = Ic*If*mVals*c10(:)./G.cells.diameters.^2;
        
        nn = size(x,1);
        mVals = bsxfun(@power, repmat([x(:,1); xq1(:,1); xq2(:,1)],1,9*9), alphaBi)...
              .*bsxfun(@power, repmat([x(:,2); xq1(:,2); xq2(:,2)],1,9*9), betaBi);
        mVals = bsxfun(@times, (mVals(ii,:) + mVals(jj,:))/12 + (mVals(nn + 1:2*nn,:) + mVals(2*nn+1:3*nn,:))*5/12, enx);
        mVals = sparseBlockDiag(mVals, nfn(f), 1);
        
        c6 = c6';
        m6cInt = Ic*If*mVals*c6(:)./G.cells.diameters.^2;
        
        c7 = c7';
        m7cInt = Ic*If*mVals*c7(:)./G.cells.diameters.^2;
        
        c9 = c9';
        m9cInt = Ic*If*mVals*c9(:)./G.cells.diameters.^2;
        
        m5cI = []; m8cI = []; m10cI = [];
        m6cI = []; m7cI = []; m9cI = [];
        
        for P = 1:G.cells.num
            
            mm = @(X) ((X(:,1)-G.cells.centroids(P,1))/G.cells.diameters(P)).^2;
            m5cI = [m5cI; polyhedronInt(G, P, mm, 2)];
            
            mm = @(X) ((X(:,2)-G.cells.centroids(P,2))/G.cells.diameters(P)).^2;
            m8cI = [m8cI; polyhedronInt(G, P, mm, 2)];
            
            mm = @(X) ((X(:,3)-G.cells.centroids(P,3))/G.cells.diameters(P)).^2;
            m10cI = [m10cI; polyhedronInt(G, P, mm, 2)];  
            
            mm = @(X) (X(:,1)-G.cells.centroids(P,1)).*(X(:,2)-G.cells.centroids(P,2))...
                      ./G.cells.diameters(P)^2;
            m6cI = [m6cI; polyhedronInt(G, P, mm, 2)];
            
            mm = @(X) (X(:,1)-G.cells.centroids(P,1)).*(X(:,3)-G.cells.centroids(P,3))...
                      ./G.cells.diameters(P)^2;
            m7cI = [m7cI; polyhedronInt(G, P, mm, 2)];
            
            mm = @(X) (X(:,2)-G.cells.centroids(P,2)).*(X(:,3)-G.cells.centroids(P,3))...
                      ./G.cells.diameters(P)^2;
            m9cI = [m9cI; polyhedronInt(G, P, mm, 2)];
%             
        end 
        
        norm(m5cI - m5cInt)
        norm(m6cI - m6cInt)
        norm(m7cI - m7cInt)
        norm(m8cI - m8cInt)
        norm(m9cI - m9cInt)
        norm(m10cI - m10cInt)
        a = 1;
        
        
        alpha = [0 1 0 0 2 1 1 0 0 0];
        beta  = [0 0 1 0 0 1 0 2 1 0];
        gamma = [0 0 0 1 0 0 1 0 1 2];
        
        ncn = diff(G.cells.nodePos);
        nce = diff(G.cells.edgePos);
        ncf = diff(G.cells.facePos);
        
        
        xMon = bsxfun(@rdivide, G.nodes.coords(G.cells.nodes,:) ...
                              - rldecode(G.cells.centroids, ncn), ...
                                rldecode(G.cells.diameters, ncn, 1));
                            
        mn =   bsxfun(@power, repmat(xMon(:,1),1, nk), alpha) ...
             .*bsxfun(@power, repmat(xMon(:,2),1, nk), beta ) ...
             .*bsxfun(@power, repmat(xMon(:,3),1, nk), gamma );
         
        ecMon = bsxfun(@rdivide, G.edges.centroids(G.cells.edges,:) ...
                              - rldecode(G.cells.centroids, nce), ...
                                rldecode(G.cells.diameters, nce, 1));
        
        me =   bsxfun(@power, repmat(ecMon(:,1),1, nk), alpha) ...
             .*bsxfun(@power, repmat(ecMon(:,2),1, nk), beta ) ...
             .*bsxfun(@power, repmat(ecMon(:,3),1, nk), gamma );
        
        fcMon = bsxfun(@rdivide, G.faces.centroids(f,:) ...
                              - rldecode(G.cells.centroids, ncf), ...
                                rldecode(G.cells.diameters, ncf, 1));
         
        m2m4fInt = bsxfun(@power, repmat(fcMon(:,1),1,3), alpha(2:4) ) ...
                 .*bsxfun(@power, repmat(fcMon(:,2),1, 3), beta (2:4) ) ...
                 .*bsxfun(@power, repmat(fcMon(:,3),1, 3), gamma(2:4) );

        dof = [0; cumsum(NP(1:end-1))] + 1;
%         sub = rldecode(-1:-1:-G.cells.num);
             
        iiN = mcolon(dof, dof + ncn - 1);
        iiE = mcolon(dof + ncn, dof + ncn + nce - 1);
        iiF = mcolon(dof + ncn + nce, dof + ncn + nce + ncf - 1);
        iiP = mcolon(dof + ncn + nce + ncf, dof + ncn + nce + ncf);
        
        D([iiN, iiE, iiF, iiP], :) ...
            = [mn; me; ones(numel(f), 1), m2m4fInt, ...
               bsxfun(@rdivide, [m5fInt, m6fInt, m7fInt, m8fInt, m9fInt, m10fInt], G.faces.areas(f)); ...
               ones(G.cells.num,1), zeros(G.cells.num, 3), ...
               bsxfun(@rdivide, [m5cInt, m6cInt, m7cInt, m8cInt, m9cInt, m10cInt], G.cells.volumes)];
        D = sparseBlockDiag(D, NP, 1);

    end
    
    %% CALCULATE B MATRICES
    
%     a5 = 1;     b5 = 0;     g5 = 0;     c5 = 2*G.faces.normals(f,1);
%     a6 = [0 1]; b6 = [1 0]; g6 = [0 0]; c6 = G.faces.normals(f,1:2);
%     a7 = [0 1]; b7 = [0 0]; g7 = [1 0]; c7 = G.faces.normals(f,[1,3]);
%     a8 = 0;     b8 = 1;     g8 = 0;     c8 = 2*G.faces.normals(f,2);
%     a9 = [0 0]; b9 = [0 1]; g9 = [1 0]; c9 = G.faces.normals(f,2:3);
%     a10 = 0; b10 = 0; g10 = 1; c10 = 2*G.faces.normals(f,3);
    
    NF = diff(G.faces.nodePos) + diff(G.faces.edgePos)*polyDim(k-2, 1) + polyDim(k-2, 2);
    
    alpha = [0 1 0 0 2 1 1 0 0 0];
    beta  = [0 0 1 0 0 1 0 2 1 0];
    gamma = [0 0 0 1 0 0 1 0 1 2];
    c     = [1 1 1 1 1 1 1 1 1 1];    
    [c_1, c_2, c_3, alphax, betay, gammaz] = polyGrad(c, alpha, beta, gamma);
    
    Kmat  = reshape(K', 3, [])';
    
    cc = Kmat*[c_1; c_2; c_3];
    
    vec = rldecode((1:3:3*G.cells.num)', diff(G.cells.facePos), 1);
    ii = mcolon(vec, vec + 2);    
    
    fn = bsxfun(@rdivide, G.faces.normals(f,:),G.faces.areas(f));
    fSign = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
                  ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
    fn = bsxfun(@times, fn, fSign);
    
    fn = reshape(fn', [], 1);
    
    cc = bsxfun(@times, cc(ii,:), fn);
    cc = reshape(cc', 3*10, [])';
    
     alpha = [alphax, alpha, alpha ];
     beta  = [beta  , betay, beta  ];
     gamma = [gamma , gamma, gammaz];
    
%     ii = ~all(cc==0, 1);
%     cc = cc(:, ii);
%     alpha = alpha(ii);
%     beta  = beta(ii);
%     gamma = gamma(ii);
    
%     cx = cx(:,[3,1,2]);
%     cy = cy(:,[3,1,2]);
%     cz = cz(:,[3,1,2]);
    
    cx = trinomialExpansion(v1(1,f)', v2(1,f)', G.faces.centroids(f,1)-ccf(:,1), 1);
    cy = trinomialExpansion(v1(2,f)', v2(2,f)', G.faces.centroids(f,2)-ccf(:,2), 1);
    cz = trinomialExpansion(v1(3,f)', v2(3,f)', G.faces.centroids(f,3)-ccf(:,3), 1);    
    
    c2 = expandGradients(alpha, beta, gamma, cc, cx, cy, cz, nk, 2)
    c3 = expandGradients(alpha, beta, gamma, cc, cx, cy, cz, nk, 3)
    c4 = expandGradients(alpha, beta, gamma, cc, cx, cy, cz, nk, 4)
    c5 = expandGradients(alpha, beta, gamma, cc, cx, cy, cz, nk, 5)
%     
    cx = cx(:, [3,1,2]); cy = cy(:,[3,1,2]); cz = cz(:,[3,1,2]);
    
    c2 = [cc(:,2), zeros(numel(f),5)];
    c3 = [cc(:,13), zeros(numel(f),5)];
    c4 = [cc(:,24), zeros(numel(f),5)];
    
%     c5 = trinomialExpansion(v1(1,f)',v2(1,f)', G.faces.centroids(f,1) - ccf(:,1), 1);
%     cy = trinomialExpansion(v1(2,f)',v2(2,f)', G.faces.centroids(f,2) - ccf(:,2), 1);
%     cz = trinomialExpansion(v1(3,f)',v2(3,f)', G.faces.centroids(f,3) - ccf(:,3), 1);
    c5 = [bsxfun(@times, cc(:,5), cx), zeros(numel(f), 3)];
    c6 = [bsxfun(@times, cc(:,6), cx) + bsxfun(@times, cc(:,16), cy), zeros(numel(f), 3)];
    c7 = [bsxfun(@times, cc(:,7), cx) + bsxfun(@times, cc(:,27), cz), zeros(numel(f), 3)];    
    c8 = [bsxfun(@times, cc(:,18), cy), zeros(numel(f), 3)];
    c9 = [bsxfun(@times, cc(:,19), cy) + bsxfun(@times, cc(:,29), cz), zeros(numel(f), 3)];
    c10 = [bsxfun(@times, cc(:,30), cz), zeros(numel(f), 3)];
%     
    
    
      
    PiNFstar = squeezeBlockDiag(PiNFstar', NF(f), sum(NF(f)), polyDim(2,2));
% 
%     
%     cx = bsxfun(@times, [cx(:,3), cx(:, 1:2), zeros(numel(f), 3)], 2*G.faces.normals(f,1));
    
%     
    a = [0 1 0 2 1 0]; b = [0 0 1 0 1 2];
    
%     c2 = [c2(:, 1:3) + c2(:, 4:6) + c2(:, 7:9), zeros(numel(f), 
    c2 = rldecode(c2, NF(f), 1);
    [alpha, beta, c2] = polyProducts(c2, PiNFstar, a, b);
    
%     c3 = [c3(:, 1:3) + c3(:, 4:6) + c3(:, 7:9), zeros(numel(f), 3)];
    c3 = rldecode(c3, NF(f), 1);
    [~, ~, c3] = polyProducts(c3, PiNFstar, a, b);
    
%     c4 = [c4(:, 1:3) + c4(:, 4:6) + c4(:, 7:9), zeros(numel(f), 3)];
    c4 = rldecode(c4, NF(f), 1);
    [~, ~, c4] = polyProducts(c4, PiNFstar, a, b);
    
%     c5 = [c5(:, 1:3) + c5(:, 4:6) + c5(:, 7:9), zeros(numel(f), 3)];
    c5 = rldecode(c5, NF(f), 1);
    [~, ~, c5] = polyProducts(c5, PiNFstar, a, b);
    
%     c6 = [c6(:, 1:3) + c6(:, 4:6) + c6(:, 7:9), zeros(numel(f), 3)];
    c6 = rldecode(c6, NF(f), 1);
    [~,~, c6] = polyProducts(c6, PiNFstar, a, b);
    
%     c7 = [c7(:, 1:3) + c7(:, 4:6) + c7(:, 7:9), zeros(numel(f), 3)];
    c7 = rldecode(c7, NF(f), 1);
    [~,~, c7] = polyProducts(c7, PiNFstar, a,b );
    
%     c8 = [c8(:, 1:3) + c8(:, 4:6) + c8(:, 7:9), zeros(numel(f), 3)];
    c8 = rldecode(c8, NF(f), 1);
    [~,  ~, c8] = polyProducts(c8, PiNFstar, a,b);
    
%     c9 = [c9(:, 1:3) + c9(:, 4:6) + c9(:, 7:9), zeros(numel(f), 3)];
    c9 = rldecode(c9, NF(f), 1);
    [~,~, c9] = polyProducts(c9, PiNFstar, a,b);
    
%     c10 = [c10(:, 1:3) + c10(:, 4:6) + c10(:, 7:9), zeros(numel(f), 3)];
    c10 = rldecode(c10, NF(f), 1);
    [~,~, c10] = polyProducts(c10, PiNFstar, a,b);
    
    
    alpha = alpha + 1;
    c2 = bsxfun(@rdivide, c2, alpha);
    c3 = bsxfun(@rdivide, c3, alpha);
    c4 = bsxfun(@rdivide, c4, alpha);
    c5 = bsxfun(@rdivide, c5, alpha);
    c6 = bsxfun(@rdivide, c6, alpha);
    c7 = bsxfun(@rdivide, c7, alpha);
    c8 = bsxfun(@rdivide, c8, alpha);
    c9 = bsxfun(@rdivide, c9, alpha);
    c10 = bsxfun(@rdivide, c10, alpha);
    
%     
%     
%     cx = bsxfun(@times, [cx(:,3), cx(:, 1:2), zeros(numel(f), 3)], 2*G.faces.normals(f,1));
%     cx = rldecode(cx, NF(f), 1);
%     
%     alpha = [0 1 0 2 1 0]; beta = [0 0 1 0 1 2];
%     [alpha, beta, c] = polyProducts(cx, PiNFstar, alpha, beta);
    
    vec = [0; cumsum(nfe(f))];
    iifn = mcolon(rldecode(vec(1:end-1)+1, NF(f), 1), rldecode(vec(2:end), NF(f), 1));
    
    pos = [1;cumsum(nfn(f))+1];
    ii = 1:size(x,1); jj = ii;
    jj(1:end-1) = jj(2:end);
    jj(cumsum(pos(2:end)-pos(1:end-1))) = ii([1;cumsum(pos(2:end-1) - pos(1:end-2))+1]);
    
    
    cDiam = rldecode(rldecode(G.cells.diameters, ncf, 1), nfn(f), 1);
    x = bsxfun(@rdivide, x, cDiam);
    ec = bsxfun(@rdivide, ec, cDiam);
    
    mVals = bsxfun(@power, repmat([x(:,1); ec(:,1)],1,6*6), alpha)...
          .*bsxfun(@power, repmat([x(:,2); ec(:,2)],1,6*6), beta);
      
%     mVals = bsxfun(@times, 
    mVals = bsxfun(@times, (mVals(ii,:) + mVals(jj,:))/6 + mVals(nn + 1:end,:)*2/3, enx);
    
    mVals = mVals(iifn,:);
    
    
    mVals = sparseBlockDiag(mVals, rldecode(nfe(f), NF(f),1), 1);
    

    
    NFf = NF(f);
    dof = [0; cumsum(NFf(1:end-1))] + 1;
    
    iiN = mcolon(dof, dof + nfn(f) - 1);
    iiE = mcolon(dof + nfn(f), dof + nfn(f) + nfe(f) - 1);
    iiF = mcolon(dof + nfn(f) + nfe(f), dof + nfn(f) + nfe(f));
    
    fDof([iiN, iiE, iiF]) = [G.faces.nodes(n,1); ...
                                  G.faces.edges(e,1) + G.nodes.num; ...
                                   f + G.nodes.num + G.edges.num];
%                                   rldecode((1:G.cells.num)', diff(G.cells.facePos),1) ...
%                                   + G.nodes.num + G.edges.num + G.faces.num]; 
    
    
    
    f2glob = sparse(fDof,1:numel(fDof),1);
    
    clear dofVec
    
    dof = [0; cumsum(NP(1:end-1))]+1;
    iiN = mcolon(dof, dof + ncn -1);
    iiE = mcolon(dof + ncn, dof + ncn + nce -1);
    iiF = mcolon(dof + ncn + nce, dof + ncn + nce + ncf - 1);
    
    dofVec([iiN, iiE, iiF]) = [G.cells.nodes; ...
                                    G.cells.edges + G.nodes.num; ...
                                    G.cells.faces(:,1) + G.nodes.num + G.edges.num];
    dofVec = dofVec(dofVec~=0);
    glob2c = sparse(1:numel(dofVec), dofVec, 1);
    
    f2c = glob2c*f2glob;
    
    If = sparseBlockDiag(ones(1, sum(rldecode(nfe(f), NF(f),1))), rldecode(nfe(f), NF(f),1), 2); 
    
    c2 = c2';
    m2Int = f2c*If*mVals*c2(:);
    
    c3 = c3';
    m3Int = f2c*If*mVals*c3(:);
        
    c4 = c4';
    m4Int = f2c*If*mVals*c4(:);
        
    c5 = c5';
    m5Int = f2c*If*mVals*c5(:);
    
    c6 = c6';
    m6Int = f2c*If*mVals*c6(:);
        
    c7 = c7';
    m7Int = f2c*If*mVals*c7(:);
        
    c8 = c8';
    m8Int = f2c*If*mVals*c8(:);
    
    c9 = c9';
    m9Int = f2c*If*mVals*c9(:);
        
    c10 = c10';
    m10Int = f2c*If*mVals*c10(:);
    
    BT = zeros(sum(NP), nk);
    
    vec = [0; cumsum(NP)] + 1;
    ii = mcolon(vec(1:end-1), vec(2:end)-2);
    
    BT(ii,2:end) = [m2Int, m3Int, m4Int, m5Int, m6Int, m7Int, m8Int, m9Int, m10Int];
    
    vec = zeros(G.cells.num,nk);
    vec(:, [1,5:nk]) = [ones(G.cells.num,1), ...
                       bsxfun(@times, -2*[K(:,1:3), K(:,5:6), K(:,9)], ...
                       G.cells.volumes./G.cells.diameters.^2)];
    BT(cumsum(NP),:) = BT(cumsum(NP),:) + vec;
    
    B = sparseBlockDiag(BT', NP, 2);
    
    M = B*D;
    
    PiNstar = M\B;
    PiN = D*PiNstar;
    
    end
    
    %%  CALCULATE B MATRICES
%     
%     if k == 1
%         Kmat = permTensor(K, 3);
%     end
    
end

    
%--------------------------------------------------------------------------

function [B, D, B1, D1] = computeBD2D(cc, cd, cv, ncn, ncf, fn, fa, fPos, x, nn, nPos, K, NP, k, dim)

B1 = []; D1 = [];

[m, grad_m, int_m] = retrieveMonomials(2,k);

nk = (k+1)*(k+2)/2;
nkk = k*(k+1)/2;

nc = numel(cv);
nf = numel(fa);
% K = zeros(2*G.cells.num, 2);
% K([1:2:end,2:2:end],:) = [K(:,1:2);K(:,3:4)];

%%  CALCULATE D MATRICES

if k == 1
    
    xMon = bsxfun(@rdivide, ...
                  x ...
                  - rldecode(cc, ncn, 1), ...
                    rldecode(cd, ncn,1));
    D = m(xMon);
    
    if dim == 2 && k == 1
        D1 = D(:, 1:nkk);
        D1 = sparseBlockDiag(D1, NP, 1);
    end
    D = sparseBlockDiag(D, ncn, 1);
    

else
    
    xMon = bsxfun(@rdivide, ...
                  x ...
                   - repmat(rldecode(cc, ncn, 1),2,1), ...
                     repmat(rldecode(cd, ncn,1),2,1));
    ii = 1:nn; jj = ii;
    jj(1:end-1) = jj(2:end);
    jj(cumsum(ncn)) = ii([1;cumsum(ncn(1:end-1))+1]);
    intD = bsxfun(@times, ( int_m(xMon(ii,:)) + int_m(xMon(jj,:)) )/6 ...
                          + int_m(xMon(nn+1:end,:))*2/3, ...
                          fn(:,1)                   );
    I = sparseBlockDiag(ones(1, sum(ncf)), ncf, 2);
    intD = bsxfun(@times, I*intD, cd./cv);
    
    
    nodeDof = mcolon([1;cumsum(NP(1:end-1))+1],[1;cumsum(NP(1:end-1))+1]+ncn-1);
    edgeDof = nodeDof + rldecode(ncn, ncn, 1)';
    D = zeros(sum(NP), nk);
    D([nodeDof, edgeDof],:) = m(xMon);
    D(cumsum(NP),:) = intD;
    
    D = sparseBlockDiag(D, NP, 1);
    
end

%%  CALCULATE B MATRICES

if k == 1
    
%     K = reshape(K, dim, [])';
    K = sparseBlockDiag(K, 2*ones(1,nc), 1);
    
    fn = bsxfun(@rdivide, fn, ...
                                 rldecode(cd, ncf, 1));
    fn = sparseBlockDiag(fn, ncf, 1);

    B = .5*fn*K;

    ii = 1:nf; jj = ii;
    jj(2:end) = jj(1:end-1);
    jj([1;cumsum(fPos(2:end-1)    ...
                -fPos(1:end-2))+1]) ...
     = ii(cumsum(fPos(2:end)-fPos(1:end-1)));

    B = B(ii,:) + B(jj,:);
    B = [.5*(fa(ii) + fa(jj)), ...
                                   squeezeBlockDiag(B,ncf, nf, nk-1)];
    
    if dim == 2 && k == 1
        B1 = B(:,1:nkk);
        B1 = sparseBlockDiag(B1', NP, 2);
    end
    
    B = sparseBlockDiag(B', NP, 2);
    
else
    
    gm = grad_m(xMon);
    ii = repmat((1:5*sum(ncn+ ncf))',2,1);
    jj = repmat(1:2,5*sum(ncn + ncf),1);
    add = repmat([rldecode((0:2:2*(nc-1))', ncn,1);rldecode((0:2:2*(nc-1))', ncf,1)], 5,1);
    jj = bsxfun(@plus,jj,add);
    jj = jj(:);
    
    intB = sparse(ii, jj, gm(:), 5*sum(ncn+ ncf), 2*nc)*K;
    
    %   Dot product by length-weighted face normals.
    
    ii = 1:nn; jj = ii;
    jj(2:end) = jj(1:end-1);
    jj(nPos(1:end-1)) = ii(nPos(2:end)-1);
    
    iin = repmat(1:nn, 1, 5) + rldecode(0:(nn+nf):4*(nn+nf), nn*ones(1,5), 2);
    iif =  iin + nn;
    
    intB = sum(intB([iin, iin, iif],:).*...
               [repmat(fn,5,1); repmat(fn(jj,:),5,1); repmat(fn,5,1)], 2);
    
    %   Evaluate line integrals using three-point Gauss-Lobatto.
           
    intB = [reshape((intB(1:numel(iin)) + intB(numel(iin)+1:2*numel(iin)))/6, nn, 5);
            reshape(intB(2*numel(iin)+1:end)*2/3, nn, 5)];
    
    intB = bsxfun(@rdivide, intB, repmat(rldecode(cd,diff(nPos),1),2,1));

    %   Assmble matrices.
    
    B = zeros(sum(NP), nk);
    B([nodeDof, edgeDof],2:nk) = intB;
    
    K = reshape(K', 4, [])';
    
    vec = zeros(nc,6);
    vec(:, [1,4:6]) = [ones(nc,1), ...
                       bsxfun(@times, -2*[K(:,1),K(:,2), K(:,4)], ...
                       cv./cd.^2)];
    B(cumsum(NP),:) = vec;

    B = sparseBlockDiag(B', NP, 2);

end

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



function [alpha, beta, coeff] = polyProducts(coeff1,coeff2,alph, bet)
    [r,c] = size(coeff1);
    cPos  = 1:c:c*c+1;
    coeff = zeros(r, cPos(end)-1);
    alpha = zeros(1, cPos(end)-1);
    beta  = zeros(1, cPos(end)-1);
    for i = 1:c
        coeff(:, cPos(i):cPos(i+1)-1) = coeff1(:, [i:end, 1:i-1]).*coeff2;
        alpha(cPos(i):cPos(i+1)-1) = alph + alph([i:end, 1:i-1]);
        beta(cPos(i):cPos(i+1)-1) = bet + bet([i:end, 1:i-1]);
    end
end


function [coeff1, coeff2, coeff3, alpha, beta, gamma] = polyGrad(coeff, alpha, beta, gamma)
    
    coeff3 = [];

    coeff1 = bsxfun(@times, coeff, alpha);
    alpha(alpha > 0) = alpha(alpha > 0) - 1;
    coeff2 = bsxfun(@times, coeff, beta);
    beta(beta > 0) = beta(beta > 0) -1;
    if ~isempty(gamma)
        coeff3 = bsxfun(@times, coeff, gamma);
        gamma(gamma > 0) = gamma(gamma > 0) -1;
    end

end

function [c, alpha] = polyAntiDerivative(c, alpha, dim)
    alpha = alpha +1;
    c = bsxfun(@rdivide, c, alpha);
end


function [c, alpha, beta] = expandGradients(alpha, beta, gamma, cc, cx, cy, cz, nk, m)
    c = [];
    alph = []; bet = [];
    
    for i = m:nk:3*nk
        
        if alpha(i) == 0 && beta(i) == 0 && gamma(i) == 0
            c = [c, cc(:, i)];
        elseif alpha(i) == 1 && beta(i) == 0 && gamma(i) == 0
            c = [c, bsxfun(@times, cc(:,i), cx)];
        elseif alpha(i) == 0 && beta(i) == 1 && gamma(i) == 0
            c = [c, bsxfun(@times, cc(:,i), cy)];
        elseif alpha(i) == 0 && beta(i) == 0 && gamma(i) == 1
            c = [c, bsxfun(@times, cc(:,i), cz)];
        else
            c = [c, cc(:,i)]
        end
        
        if alpha(i) == 0 && beta(i) == 0 && gamma(i) == 0
            alph = [alph, 0]; bet = [bet, 0];
        elseif alpha(i) < 2 && beta(i) < 2 && gamma(i) < 2
            alph = [alph, 1 0 0]; bet = [bet, 0 1 0];
        else
            alph = [alph, 0]; bet = [bet, 0];
        end
        
    end
%     
    cc = zeros(size(cc,1), 3);
    cc(:,1) = sum(c(:,alph==0 & bet==0),2);
    cc(:,2) = sum(c(:,alph==1 & bet==0),2);
    cc(:,3) = sum(c(:,alph==0 & bet==1),2);
    alpha = [0 1 0 2 1 0]; beta = [0 0 1 0 1 2];
    c = cc;
end


function nk = polyDim(k, dim)
    nk = nchoosek(k+dim,k);
end

%         xExp = alpha(:,1:2:end) ~= 0;
%         yExp = ~xExp & alpha(:,2:2:end) ~=0;
%         alpha(:,[1:2:end,2:2:end]) = alpha(:, [1:2:end,2:2:end]) + [xExp, yExp];

%         xyExp = zeros(6, size(T,2));
%         xyExp(:,[1:2:end,2:2:end]) = [xExp,yExp];
%         xyExp = rldecode(reshape(xyExp',2,[])', repmat(nfn,6,1),1);