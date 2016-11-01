function [T, flux] = computeNonLinearTrans(G, coSet, cellPressure) 
    if min(double(cellPressure)) < 0
        warning('Negative pressure in cells. Will fall back to linear TPFA.');
    end
    flux = nan;
    T = computeJumpTransmissibilities(G, coSet, cellPressure);
    T_face = computeContTransmissibilities(G, coSet, cellPressure);
    
    jump = coSet.jumpFace;
    for i = 1:2
        T{i} = T{i}.*jump + T_face{i}.*~jump;
    end
end

function d = computeWeight(lnorm, c, p, exclude)
    d = 0;
    for i = 1:numel(p)
        d = d + p{i}.*c(:, i).*(~exclude(:, i));
    end
    d = d.*lnorm;
end

function T = computeJumpTransmissibilities(G, coSet, cellPressure)
    cellNo = getCellNoFaces(G);
    faceNo = G.cells.faces(:, 1);
    
    dim = G.griddim;
    n_hf = size(faceNo, 1);
    lnorm = sqrt(sum(coSet.l.^2, 2));

    
    pp = cell(dim, 1);
    [pp{:}] = deal(double2ADI(zeros(n_hf, 1), cellPressure));
    [p_c, p_f] = deal(pp);
    for i = 1:dim
        p_c{i} = coSet.pressureOperators{i, 1}(cellPressure);
        p_f{i} = coSet.pressureOperators{i, 2}(cellPressure);
    end
    
    Ac = coSet.active{1};
    Af = coSet.active{2};
    d_c = computeWeight(lnorm, coSet.C{1}, p_c, Ac);
    d_f = computeWeight(lnorm, coSet.C{2}, p_f, Af);
    
    thold = 0;
    bad = d_c <= thold & d_f <= thold;
    d_c(bad) = 1;
    d_f(bad) = 1;
    
    d_tot = d_c + d_f;
    
    % intentional swapping of terms
    u_c = d_f./d_tot;
    u_f = d_c./d_tot;
    
    deviation = -u_f.*d_f + u_c.*d_c;
    unity = u_f + u_c;
    
    all_c = coSet.C{1};
    all_f = coSet.C{2};
    
%     [max(abs(deviation.val)), max(abs(double(unity)- 1))]
    N_c = lnorm.*(u_c.*(sum(all_c, 2)) + u_f.*(sum(Af.*all_f, 2)));
    N_f = lnorm.*(u_f.*(sum(all_f, 2)) + u_c.*(sum(Ac.*all_c, 2)));

    H = sparse(faceNo, (1:numel(faceNo))', 1);
    N_tot = H*N_f;
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

        for i = 1:2
            T{i} = T{i}./N_tot;
        end
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


function [T, intx] = computeContTransmissibilities(G, coSet, cellPressure)
    intx = find(all(G.faces.neighbors > 0, 2));
    
    dim = G.griddim;
    n_intf = size(intx, 1);
    lnorm = sqrt(sum(coSet.faceSet.l.^2, 2));
   
    
    pp = cell(dim, 1);
    [pp{:}] = deal(double2ADI(zeros(n_intf, 1), cellPressure));
    [p_l, p_r] = deal(pp);
    for i = 1:dim
        p_l{i} = coSet.faceSet.pressureOperators{i, 1}(cellPressure);
        p_r{i} = coSet.faceSet.pressureOperators{i, 2}(cellPressure);
    end
    
    A_l = coSet.faceSet.active{1};
    A_r = coSet.faceSet.active{2};
    d_l = computeWeight(lnorm, coSet.faceSet.C{1}, p_l, A_l);
    d_r = computeWeight(lnorm, coSet.faceSet.C{2}, p_r, A_r);
    
    thold = 0;
    bad = d_l <= thold & d_r <= thold;
    d_l(bad) = 1;
    d_r(bad) = 1;
    
    d_tot = d_l + d_r;
    
    % intentional swapping of terms
    u_l = d_r./d_tot;
    u_r = d_l./d_tot;
    
    % deviation = -u_r.*d_r + u_l.*d_l;
    % unity = u_r + u_l;
    
    all_l = coSet.faceSet.C{1};
    all_r = coSet.faceSet.C{2};
    
    T_l = lnorm.*(u_l.*(sum(all_l, 2)) + u_r.*(sum(A_r.*all_r, 2)));
    T_r = lnorm.*(u_r.*(sum(all_r, 2)) + u_l.*(sum(A_l.*all_l, 2)));

    TT = (double2ADI(zeros(G.faces.num, 1), cellPressure));
    T = {TT, TT};
    T{1}(intx) = T_l;
    T{2}(intx) = T_r;
end