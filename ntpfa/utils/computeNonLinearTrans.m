function [T, flux] = computeNonLinearTrans(G, coSet, cellPressure)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if min(cellPressure.value) < 0
        warning('Negative pressure in cells. Will fall back to linear TPFA.');
    end
    flux = nan;
    T_hf = computeJumpTransmissibilities(G, coSet, cellPressure);
    T_face = computeContTransmissibilities(G, coSet, cellPressure);
    
    jump = coSet.jumpFace;
    T = T_face;
    for i = 1:2
        T{i} = T_hf{i}.*jump + ~jump.*T_face{i};
    end
end

function d = computeWeight(lnorm, c, p, exclude)
    d = 0;
    for i = 1:numel(p)
        d = d + p{i}.*c(:, i).*(~exclude(:, i));
    end
    d = d.*lnorm;
end

function d = computeWeight2(lnorm, c, p)
    d = 0;
    for i = 1:numel(p)
        d = d + p{i}.*c(:, i);
    end
    d = d.*lnorm;
end

function [u_c, u_f] = get_mu(d_c, d_f)
    thold = 0;
    bad = d_c <= thold | d_f <= thold;
    d_c(bad) = 1;
    d_f(bad) = 1;
    
    d_tot = d_c + d_f;
    % intentional swapping of terms
    u_c = d_f./d_tot;
    u_f = d_c./d_tot;
end

function T = computeJumpTransmissibilities(G, coSet, cellPressure)
    cellNo = getCellNoFaces(G);
    faceNo = G.cells.faces(:, 1);
    
    dim = G.griddim;
    n_hf = size(faceNo, 1);
    lnorm = sqrt(sum(coSet.l.^2, 2));
    
    elimWeights = false;
    
    pp = cell(dim, 1);
    [pp{:}] = deal(double2ADI(zeros(n_hf, 1), cellPressure));
    [p_c, p_f] = deal(pp);
    
    sum_c = zeros(n_hf, G.griddim);
    sum_f = zeros(n_hf, G.griddim);
    for i = 1:dim
        op_c = coSet.pressureOperators{i, 1};
        op_f = coSet.pressureOperators{i, 2};
        
        if elimWeights
            p_c{i} = op_c.P_passive*cellPressure;
            p_f{i} = op_f.P_passive*cellPressure;
        
            sum_c(:, i) = sum(op_c.P_active, 2);
            sum_f(:, i) = sum(op_f.P_active, 2);
        else
            p_c{i} = op_c.P*cellPressure;
            p_f{i} = op_f.P*cellPressure;
        end
    end
    
    Ac = coSet.active{1};
    Af = coSet.active{2};
    if elimWeights
        d_c = computeWeight2(lnorm, coSet.C{1}, p_c);
        d_f = computeWeight2(lnorm, coSet.C{2}, p_f);
    else
        d_c = computeWeight(lnorm, coSet.C{1}, p_c, Ac);
        d_f = computeWeight(lnorm, coSet.C{2}, p_f, Af);
    end

    [u_c, u_f] = get_mu(d_c, d_f);

    deviation = -u_f.*d_f + u_c.*d_c;
    unity = u_f + u_c;
    
    all_c = coSet.C{1};
    all_f = coSet.C{2};
    
    if elimWeights
        N_c = lnorm.*(u_c.*(sum(all_c, 2)) + u_f.*(sum(all_f.*sum_f, 2)));
        N_f = lnorm.*(u_f.*(sum(all_f, 2)) + u_c.*(sum(all_c.*sum_c, 2)));
    else
        N_c = lnorm.*(u_c.*(sum(all_c, 2)) + u_f.*(sum(all_f.*Af, 2)));
        N_f = lnorm.*(u_f.*(sum(all_f, 2)) + u_c.*(sum(all_c.*Ac, 2)));
    end
    
    H = sparse(faceNo, (1:numel(faceNo))', 1);
    N_tot = H*N_f;
    assert(all(N_tot >= 0))
    
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
    assert(all(N_f >= 0));
    assert(all(T{1}.value>=0))
    assert(all(T{2}.value>=0))
end


function [T, intx] = computeContTransmissibilities(G, coSet, cellPressure)
    intx = find(all(G.faces.neighbors > 0, 2));
    
    dim = G.griddim;
    n_intf = size(intx, 1);
    lnorm = sqrt(sum(coSet.faceSet.l.^2, 2));
   
    elimWeights = false;
    
    pp = cell(dim, 1);
    [pp{:}] = deal(double2ADI(zeros(n_intf, 1), cellPressure));
    [p_l, p_r] = deal(pp);
    sum_l = zeros(n_intf, G.griddim);
    sum_r = zeros(n_intf, G.griddim);

    for i = 1:dim
        op_l = coSet.faceSet.pressureOperators{i, 1};
        op_r = coSet.faceSet.pressureOperators{i, 2};
        
        if elimWeights
            p_l{i} = op_l.P_passive*cellPressure;
            p_r{i} = op_r.P_passive*cellPressure;
            sum_l(:, i) = sum(op_l.P_active, 2);
            sum_r(:, i) = sum(op_r.P_active, 2);
        else
            p_l{i} = op_l.P*cellPressure;
            p_r{i} = op_r.P*cellPressure;
        end
    end
    
    if elimWeights
        d_l = computeWeight2(lnorm, coSet.faceSet.C{1}, p_l);
        d_r = computeWeight2(lnorm, coSet.faceSet.C{2}, p_r);
    else
        A_l = coSet.faceSet.active{1};
        A_r = coSet.faceSet.active{2};
        d_l = computeWeight(lnorm, coSet.faceSet.C{1}, p_l, A_l);
        d_r = computeWeight(lnorm, coSet.faceSet.C{2}, p_r, A_r);
    end
    [u_l, u_r] = get_mu(d_l, d_r);
    % deviation = -u_r.*d_r + u_l.*d_l;
    % unity = u_r + u_l;
    
    all_l = coSet.faceSet.C{1};
    all_r = coSet.faceSet.C{2};
    if elimWeights
        T_l = lnorm.*(u_l.*(sum(all_l, 2)) + u_r.*(sum(sum_r.*all_r, 2)));
        T_r = lnorm.*(u_r.*(sum(all_r, 2)) + u_l.*(sum(sum_l.*all_l, 2)));
    else
        T_l = lnorm.*(u_l.*(sum(all_l, 2)) + u_r.*(sum(all_r.*A_r, 2)));
        T_r = lnorm.*(u_r.*(sum(all_r, 2)) + u_l.*(sum(all_l.*A_l, 2)));
    end
    
    TT = (double2ADI(zeros(G.faces.num, 1), cellPressure));
    T = {TT, TT};
    T{1}(intx) = T_l;
    T{2}(intx) = T_r;
end
