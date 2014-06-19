function [B, varargout] = createMSFVBasis(A, DG, useCorrection)

    n_i = length(DG.ii);
    n_e = n_i + length(DG.ee);
    n_l = n_e + length(DG.ll);
    n_n = n_l + length(DG.nn);

    if isempty(DG.ll)
        dim = 2;
    else
        dim = 3;
    end

    ind.i = 1:n_i;
    ind.e = (n_i+1):n_e;
    ind.l = n_e+1:n_l;
    ind.n = n_l+1:n_n;

    A = DG.P*A*DG.P';

    A_ii = A(ind.i, ind.i);
    A_ie = A(ind.i, ind.e);
    A_ee = A(ind.e, ind.e);

    A_ll = A(ind.l, ind.l);
    A_le = A(ind.l, ind.e);
    A_li = A(ind.l, ind.i);
    A_ln = A(ind.l, ind.n);
    A_el = A(ind.e, ind.l);


    M_ee = A_ee + diag(sum(A_ie, 1));

    if dim == 3
        M_ll = A_ll + diag(sum(A_le, 2)) + diag(sum(A_li, 2));
        B = createBasis3D(A_ie, A_ii, M_ee, M_ll, A_ln, A_el, @mldivide_scalable);
    else
        A_en = A(ind.e, ind.n);
        B = createBasis2D(A_ie, A_ii, M_ee, A_en);
    end
    B = DG.P'*B;

    if useCorrection
        invii_ie = mldivide(A_ii, A_ie);
        if dim == 3
            invee_el = mldivide_scalable(M_ee, A_el);
            Cx = @(r) DG.P'*correctionOperator3D(M_ee, A_ii, M_ll, invee_el, invii_ie, ind, DG.P*r);
        else
            Cx = @(r) DG.P'*correctionOperator2D(M_ee, A_ii, invii_ie, ind, DG.P*r);
        end
        varargout{1} = Cx;
    else
        varargout{1} = Cxr = @(r) 0*r;
    end
end

% Basis functions

function B = createBasis3D(A_ie, A_ii, M_ee, M_ll, A_ln, A_el, lsolve)
    ll_ln = (lsolve(M_ll,A_ln));
    tmp = (lsolve(M_ee,A_el*ll_ln));
    B = [-(lsolve(A_ii,(A_ie*tmp)));...
        tmp;...
        -ll_ln;...
        eye(size(ll_ln,2))];
    return
end

function B = createBasis2D(A_ie, A_ii, M_ee, A_en)
    tmp = (M_ee\A_en);
    B = [(A_ii\(A_ie*tmp));...
        -tmp;...
        eye(size(tmp,2))];
    return
end

% Correction functions

function res = correctionOperator3D(M_ee, A_ii, M_ll, invee_el, invii_ie, ind, r)
    res = zeros(size(r));
    invee_r = mldivide(M_ee, r(ind.e));
    invll_r = mldivide(M_ll, r(ind.l));

    res(ind.i) = mldivide(A_ii, r(ind.i)) - invii_ie*mldivide(M_ee, r(ind.e)) + invii_ie*invee_el*mldivide(M_ll, r(ind.l));
    res(ind.e) = invee_r - invee_el*invll_r;
    res(ind.l) = invll_r;
end

function res = correctionOperator2D(M_ee, A_ii, invii_ie, ind, r)
    res = zeros(size(r));
    invee_r = mldivide(M_ee, r(ind.e));

    res(ind.i) = mldivide(A_ii, r(ind.i)) - invii_ie*invee_r;
    res(ind.e) = invee_r;
end

function x = mldivide_scalable(A, b)
    x = [];
    max_sol = 250;
    n_rhs = size(b, 2);
    n_sub = ceil(n_rhs/max_sol);

    for i = 1:n_sub
        tmp = A\b(:, ((i-1)*max_sol + 1):min((i*max_sol), n_rhs));
        x = horzcat(x, tmp); %#ok
    end
end