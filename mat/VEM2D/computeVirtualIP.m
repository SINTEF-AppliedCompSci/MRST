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

    
%     [ii, jj] = blockDiagIndex(repmat(3, [G.cells.num,1]));
%     Kmat = rldecode(K, diff(G.cells.facePos),1);
%     Kmat = reshape(Kmat', 3, [])';
    
    
    %   Function space dimensions.
    
    
    f = G.cells.faces(:,1);
    
%     [ii,jj] = blockDiagIndex(repmat(3, [G.cells.num,1]));
%     Kmat = K';
%     Kmat = sparse(ii,jj,Kmat(:));
    
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
    PiNFstar = sparse(ii, jj, invv(full(MF(kk)), repmat(nk, [G.cells.num, 1])))*BF;

    clear BF DF;
    
    %%  CALCULATE D MATRICES
    
    [m, grad_m, int_m] = retrieveMonomials(3, k);
    ncn = diff(G.cells.nodePos);
    nce = diff(G.cells.edgePos);
    ncf = diff(G.cells.facePos);
    
    
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
        xMon = bsxfun(@rdivide, [G.nodes.coords(G.cells.nodes   , :); ...
                                 G.edges.centroids(G.cells.edges, :)] ...
                       - rldecode(repmat(G.cells.centroids,2,1), [ncn; nce],1), ...
                         rldecode(repmat(G.cells.diameters,2,1), [ncn; nce], 1));
        DfInt = m(bsxfun(@rdivide, G.faces.centroids(G.cells.faces(:,1),:) ...
                                 - rldecode(G.cells.centroids, ncf, 1), ...
                                   rldecode(G.cells.diameters, ncf, 1)));

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
        
        for P = 1:G.cells.num
            
            ff = G.cells.faces(G.cells.facePos(P):G.cells.facePos(P+1)-1);
            
            mm = @(X) ((X(:,1)-G.cells.centroids(P,1))/G.cells.diameters(P)).^2;
            m5fI = [m5fI; polygonInt3D(G, ff, mm, 2)];
            
            mm = @(X) ((X(:,2)-G.cells.centroids(P,2))/G.cells.diameters(P)).^2;
            m8fI = [m8fI; polygonInt3D(G, ff, mm, 2)];
            
            mm = @(X) ((X(:,3)-G.cells.centroids(P,3))/G.cells.diameters(P)).^2;
            m10fI = [m10fI; polygonInt3D(G, ff, mm, 2)];  
            
            mm = @(X) (X(:,1)-G.cells.centroids(P,1)).*(X(:,2)-G.cells.centroids(P,2))...
                      ./G.cells.diameters(P)^2;
            m6fI = [m6fI; polygonInt3D(G, ff, mm, 2)];
            
            mm = @(X) (X(:,1)-G.cells.centroids(P,1)).*(X(:,3)-G.cells.centroids(P,3))...
                      ./G.cells.diameters(P)^2;
            m7fI = [m7fI; polygonInt3D(G, ff, mm, 2)];
            
            mm = @(X) (X(:,2)-G.cells.centroids(P,2)).*(X(:,3)-G.cells.centroids(P,3))...
                      ./G.cells.diameters(P)^2;
            m9fI = [m9fI; polygonInt3D(G, ff, mm, 2)];
            
        end 
        
        norm(m5fI - m5fInt)
        norm(m6fI - m6fInt)
        norm(m7fI - m7fInt)
        norm(m8fI - m8fInt)
        norm(m9fI - m9fInt)
        norm(m10fI - m10fInt)
        
        
        fSign = (-ones(numel(f),1)).^(G.faces.neighbors(f,1) ...
                  ~= rldecode((1:G.cells.num)', diff(G.cells.facePos), 1)); 
        fnx = G.faces.normals(f,1).*fSign;
        
        c5 = trinomialExpansion(v1(1,f)', v2(1,f)', G.faces.centroids(f,1) - ccf(:,1), 3);
        c8 = trinomialExpansion(v1(2,f)', v2(2,f)', G.faces.centroids(f,2) - ccf(:,2), 3);
        c10 = trinomialExpansion(v1(3,f)', v2(3,f)', G.faces.centroids(f,3) - ccf(:,3), 3);
        
        c5 = bsxfun(@times, c5, fnx/3);
        c8 = bsxfun(@times, c8, fnx/3);
        c10 = bsxfun(@times, c10, fnx/3);
        
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
        
        xq1 = bsxfun(@minus, ec, .5*sqrt(1/5)*G.edges.lengths(e));
        xq2 = bsxfun(@plus, ec,  .5*sqrt(1/5)*G.edges.lengths(e));
        
        
        nn = size(x,1);
        mVals = bsxfun(@power, repmat([x(:,1); xq1(:,1); xq2(:,1)],1,10), alphaQuad)...
              .*bsxfun(@power, repmat([x(:,2); xq1(:,2); xq2(:,2)],1,10), betaQuad);
        mVals = bsxfun(@times, (mVals(ii,:) + mVals(jj,:))/12 + (mVals(nn + 1:2*nn,:) + mVals(2*nn+1:3*nn,:))*5/12, enx);
        mVals = sparseBlockDiag(mVals, nfn(f), 1);
        
        If = sparseBlockDiag(ones(1, sum(nfn(f))), nfn(f), 2); 
        Ic = sparseBlockDiag(ones(1, sum(ncf)), ncf, 2); 
        cd = rldecode(G.cells.diameters, diff(G.cells.facePos), 1);
        
        c5 = c5';
        m5cInt = Ic*If*mVals*c5(:)./G.cells.diameters.^2;
%         
%         c8 = c8';
%         m8fInt = Ic*If*mVals*c8(:)./cd.^2;
%         
%         c10 = c10';
%         m10fInt = cI*If*mVals*c10(:)./cd.^2;
%         
        m5cI = [];
        
        for P = 1:G.cells.num
            
            mm = @(X) ((X(:,1)-G.cells.centroids(P,1))/G.cells.diameters(P)).^2;
            m5cI = [m5cI; polyhedronInt(G, P, mm, 2)];
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
        end 
        
        norm(m5cI - m5cInt)
        
        a = 1;
                
%         [alphaBi, betaBi, cBi]  = polyProducts([x; x; y], [y; z; z], G.faces.num);% 6, 7, 9
        
%         alphaBi = alphaBi+1;
%         cBi = bsxfun(@rdivide, cBi, alphaBi);
        
        
        
%         
%         monValBi = bsxfun(@power, repmat([x(:,1); ec(:,1)],1,6*6),alphaBi).*bsxfun(@power, repmat([x(:,2); ec(:,2)],1,6*6), betaBi);
%         monValBi = (bsxfun(@times, monValBi(ii,:), enx) + bsxfun(@times, monValBi(jj,:), enx))/6 ...
%                 + bsxfun(@times, monValBi(size(x,1)+1:end,:), enx)*2/3;
%         monValBi = sparseBlockDiag(monValBi, nfn(f), 1);
%         cBi      = sparseBlockDiag(cBi(f,:)', ones(1, numel(f)), 2);
%         
%         fInt = monValQuad*c5 + monValBi*cBi;
%         
%         I = sparseBlockDiag(ones(1, sum(nfn(f))), nfn(f), 2);
%         fInt = I*fInt
% 
%         cc = G.cells.centroids;
%         cn = G.cells.num;
%         cd = G.cells.diameters;
%         tic
%         m = @(x) [bsxfun(@rdivide, bsxfun(@minus, repmat(x(:,1),1,cn), cc(:,1)').^2, cd'.^2)];
% %                   bsxfun(@rdivide, bsxfun(@minus, repmat(x(:,2),1,cn), cc(:,2)').^2, cd'.^2), ...
% %                   bsxfun(@rdivide, bsxfun(@minus, repmat(x(:,3),1,cn), cc(:,3)').^2, cd'.^2)];
% %                   bsxfun(@rdivide, bsxfun(@minus, repmat(x(:,1),1,cn), cc(:,1)').* ...
% %                                    bsxfun(@minus, repmat(x(:,2),1,cn), cc(:,2)'), cd'.^2), ...
% %                   bsxfun(@rdivide, bsxfun(@minus, repmat(x(:,1),1,cn), cc(:,1)').* ...
% %                                    bsxfun(@minus, repmat(x(:,3),1,cn), cc(:,3)'), cd'.^2), ...
% 
% %                   bsxfun(@rdivide, bsxfun(@minus, repmat(x(:,2),1,cn), cc(:,2)').* ...
% %                                    bsxfun(@minus, repmat(x(:,3),1,cn), cc(:,3)'), cd'.^2), ...
% 
%               
%         fInt = polygonInt3D(G, 1:G.faces.num, m, 2);
%         cInt = polyhedronInt(G, 1:G.cells.num, m, 2);
%         toc
        
        

        
        %         

        
        

%         T = zeros(3, 2*G.faces.num);
%         T(:,[1:2:end, 2:2:end]) = [v1, v2];
%         alpha = [2 0 0; 1 1 0; 1 0 1; 0 2 0; 0 1 1; 0 0 2];
%         alpha = alpha*T;
%         
%         
%         cc = (rldecode(G.cells.centroids,ncf,1)-G.faces.centroids)*T;
% 
%         alpha(:,1:2:end) = alpha(:,1:2:end) + 1;
% %         alpha = rldecode(reshape(alpha',2,[])', repmat(nfn,6,1),1);
%         
%         alpha = reshape(alpha,6*2,[]);
%         alpha = reshape(alpha(:,f),6,[]);
%         
%         
%         
%         
%         ii = 1:sum(nfn); jj = ii;
%         jj(2:end) = jj(1:end-1);
%         jj([1;cumsum(G.faces.nodePos(2:end-1)    ...
%             -G.faces.nodePos(1:end-2))+1]) ...
%             = ii(cumsum(G.faces.nodePos(2:end)-G.faces.nodePos(1:end-1)));
%         
%         int = [repmat(x, 6, 1);repmat(ec,6,1)].^repmat(alpha,2,1);
%         int = int(:,1).*int(:,2);
%         int = (int(1:6*sum(nfn)).*repmat(en(:,1),6,1) + int(1:6*sum(nfn)).*repmat(en(jj,1),6,1))/6 ...
%             + int(6*sum(nfn)+1:end).*repmat(en(:,1),6,1)*2/3;
%         
%         I = sparseBlockDiag(ones(1, 6*sum(nfe)), repmat(nfe,6,1), 2);
%         DfInt(:,5:end) = reshape(I*int,G.faces.num, 6);
        
        
%         
%         int = 
%         int = ((int(:,1).*int(:,2)).*repmat(en(:,1),6,1) + ...
%               (int(:,1).*int(:,2)).*repmat(en(jj,1),6,1))./alpha(:,1);
% 
%         DfInt
        
%         xyInt = 
%         rldecode([xInt(:), yInt(:)], repmat(nfn,6,1),1);.*en,2)
        
%         
%         enf = sum([xInt(:),yInt(:)].*rldecode(en,2);
        
%         alpha = sparseBlockDiag(alpha*T, 2*ones(G.faces.num,1), 2);
        
        
    end
    
    %%  CALCULATE B MATRICES
%     
%     if k == 1
%         Kmat = permTensor(K, 3);
%     end
    
end

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
    
    if n == 1
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
    cPos = 1:c:c*(r+1);
    for i = 1:c
        coeff(:, cPos(i):cPos(i+1)-1) = coeff1(:, [i:end, 1:i-1]).*coeff2;
        alpha(cPos(i):cPos(i+1)-1) = alph + alph([i:end, 1:i-1]);
        beta(cPos(i):cPos(i+1)-1) = bet + bet([i:end, 1:i-1]);
    end
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