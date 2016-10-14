function [T, flux] = computeNonLinearTrans(G, coSet, cellPressure) 
    if min(double(cellPressure)) < 0
        warning('Negative pressure in cells. Will fall back to linear TPFA.');
    end
    cellNo = getCellNoFaces(G);
    faceNo = G.cells.faces(:, 1);
    
    dim = G.griddim;
    n_hf = size(faceNo, 1);
    
    intx = all(G.faces.neighbors > 0, 2);

    lnorm = sqrt(sum(coSet.l.^2, 2));
   
%     facePressure = interpolatedFacePressure(G, cellPressure, true);
%     nodePressure = nan;
    
    Tn_c = [coSet.T_norm{:, 1}];
    Tn_f = [coSet.T_norm{:, 2}];
    
    pp = cell(dim, 1);
    [pp{:}] = deal(double2ADI(zeros(n_hf, 1), cellPressure));
    [p_c, p_f] = deal(pp);
    for i = 1:dim
        p_c{i} = coSet.pressureOperators{i, 1}(cellPressure);
        p_f{i} = coSet.pressureOperators{i, 2}(cellPressure);
    end
    
    Ac = coSet.active{1};
    Af = coSet.active{2};
    d_c = computeWeight(lnorm, Tn_c, coSet.C{1}, p_c, Ac);
    d_f = computeWeight(lnorm, Tn_f, coSet.C{2}, p_f, Af);
    
%     thold = sqrt(eps);
    thold = 0;
    bad = d_c <= thold | d_f <= thold;
    d_c(bad) = 1;
    d_f(bad) = 1;
    
    d_tot = d_c + d_f;
    
    % intentional swapping of terms
    u_c = d_f./d_tot;
    u_f = d_c./d_tot;
    
    
    all_c = coSet.C{1}./Tn_c;
    all_f = coSet.C{2}./Tn_f;
    N_c = lnorm.*(u_c.*(sum(all_c, 2)) + u_f.*(sum(Af.*all_f, 2)));
    N_f = lnorm.*(u_f.*(sum(all_f, 2)) + u_c.*(sum(Ac.*all_c, 2)));

    H = sparse(faceNo, (1:numel(faceNo))', 1);
    N_tot = H*N_f;
%     N_tot = accumarray(faceNo, N_f);
    assert(all(N_tot > 0))
    
    left  = G.faces.neighbors(faceNo, 1) == cellNo;
    right = ~left;
    
    zT = double2ADI(zeros(G.faces.num, 1), cellPressure);
    T = {zT, zT};
    flux = {zT, zT};
    if 1
        T{1}(faceNo(left)) = N_c(left);
        T{1}(faceNo(right)) = T{1}(faceNo(right)).*N_f(right);

        T{2}(faceNo(right)) = N_c(right);
        T{2}(faceNo(left)) = T{2}(faceNo(left)).*N_f(left);

%         T(faceNo(right), 2) = N_c(right);
%         T(faceNo(left), 2) = T(faceNo(left), 2).*N_f(left);
        for i = 1:2
            T{i} = T{i}./N_tot;
        end
%         T = bsxfun(@rdivide, T, N_tot);
    else
        % Loop based code for debugging
        for f = 1:G.faces.num
            approxFp = zeros(G.faces.num, 1);
            subs = find(faceNo == f);
            
            if numel(subs) == 1
                continue
            end
            
            p = subs(1);
            m = subs(2);
            
            nf_p = N_f(p);
            nf_m = N_f(m);
            
            nc_p = N_c(p);
            nc_m = N_c(m);
            
            T(f, 1) = nc_p*nf_m./(nf_p + nf_m);
            T(f, 2) = nc_m*nf_p./(nf_p + nf_m);
            
            % debuggin'
            c1 = G.faces.neighbors(f, 1);
            c2 = G.faces.neighbors(f, 2);
            
            approxFp(f) = (cellPressure(c1)*nc_p + cellPressure(c2)*nc_m)./(nf_p + nf_m);
            flux(f, 1) =   nc_p*cellPressure(c1) - nf_p*approxFp(f);
            flux(f, 2) = -(nc_m*cellPressure(c2) - nf_m*approxFp(f));
            
            e = flux(f, 1) - flux(f, 2);
            if isnan(e)
                e = 0;
            end
            assert(abs(e) < 1e-12)
        end
    end
    assert(all(N_f > 0));
    assert(all(double(T{1})>=0))
    assert(all(double(T{2})>=0))
end

function d = computeWeight(lnorm, tnorm, c, p, exclude)
    d = 0;
    for i = 1:numel(p)
        d = d + p{i}.*c(:, i).*(~exclude(:, i))./tnorm(:, i);
    end
    d = d.*lnorm;
%     d = lnorm.*sum(p.*c.*(~exclude)./tnorm, 2);
end

