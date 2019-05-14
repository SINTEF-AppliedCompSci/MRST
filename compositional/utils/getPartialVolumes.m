function [V, dVdp] = getPartialVolumes(model, problem, varargin)
    state = problem.state;
    opt = struct('units', 'mass', ...
                 'pressureEpsilon', 0.1*psia, ...
                 'singlePhaseStrategy', 'numerical', ...
                 'twoPhaseStrategy', [], ...
                 'iteration', nan, ...
                 'singlePhaseDifferentiation', 'numerical', ...
                 'twoPhaseDifferentiation', []);
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.twoPhaseDifferentiation)
        opt.twoPhaseDifferentiation = opt.singlePhaseDifferentiation;
    end
    if isempty(opt.twoPhaseStrategy)
        opt.twoPhaseStrategy = opt.singlePhaseStrategy;
    end
    % Choices for strategy:
    % - Numerical (compute from matrix)
    % - Analytical (1/rho). At the moment, only for single-phase
    % - EOS (use flash)
    %
    % Differentation:
    % - Numerical
    % - Perturbation (Only valid for EOS at the moment)
    % - Analytical (Only for single-phase)
    % - None
    if nargout == 1
        opt.twoPhaseDifferentiation    = 'none';
        opt.singlePhaseDifferentiation = 'none';
    end
    
    ncell = numelValue(state.pressure);
    ncomp = model.water + numel(model.EOSModel.fluid.names);
    
    [liquid, vapor, twoPhase] = model.getFlag(state);
    singlePhase = liquid | vapor;
    
    equalStrategy = strcmpi(opt.twoPhaseDifferentiation, opt.singlePhaseDifferentiation) && ...
                    strcmpi(opt.singlePhaseStrategy, opt.twoPhaseStrategy);
    
    onlySingle = all(singlePhase);
    onlyMulti = all(twoPhase);
    uniform = onlySingle || onlyMulti || equalStrategy;
    
    if uniform
        % All cells are either in the same phase state, or we have the same
        % strategy everywhere. We can do a single call to get the correct
        % values.
        if equalStrategy || onlySingle
            s = opt.singlePhaseStrategy;
            ds = opt.singlePhaseDifferentiation;
        else
            s = opt.twoPhaseStrategy;
            ds = opt.twoPhaseDifferentiation;
        end
        subs = true(ncell, 1);
        [V, dVdp] = getPartialVolumesWrapper(model, problem, opt, subs, s, ds);
    else
        V = nan(ncell, ncomp);
        dVdp = nan(ncell, ncomp);
        % Treat single-phase
        [V(singlePhase, :), dVdp(singlePhase, :)] = getPartialVolumesWrapper(model, problem, opt, singlePhase, ...
                                    opt.singlePhaseStrategy, opt.singlePhaseDifferentiation);
        % Treat two-phase
        assert(~strcmpi(opt.twoPhaseStrategy, 'analytical'));
        [V(twoPhase, :), dVdp(twoPhase, :)] = getPartialVolumesWrapper(model, problem, opt, twoPhase, ...
                                    opt.twoPhaseStrategy, opt.twoPhaseDifferentiation);
                                
    end
end

function [V, dVdp] = getPartialVolumesWrapper(model, problem, opt, subs, strategy, dstrategy)
    state = problem.state;
    if ~all(subs)
        state = makeSubstate(state, subs);
    end
    switch lower(strategy)
        case 'numerical'
            V = getWeights(problem, subs);
        case 'eos'
            computeDerivatives = strcmpi(dstrategy, 'eos');
            [V, dVdp] = getPartialVolumesInternal(model, state, opt, subs, computeDerivatives);
            if ~computeDerivatives
                % Analytical weights not possible
                assert(strcmpi(dstrategy, 'none') || strcmpi(dstrategy, 'numerical'));
            end
        case 'analytical'
            computeDerivatives = strcmpi(dstrategy, 'analytical');
            pureLiquid = state.L > 0.5;
            pureVapor = not(pureLiquid);
            [V, dVdp] = getSinglePhaseVolumes(model, state, pureLiquid, pureVapor, computeDerivatives);
        otherwise
            error('Unsupported strategy %s.', strategy);
    end
    switch lower(dstrategy)
        case 'numerical'
            if opt.iteration > 1 && isfield(state, 'w_p')
                dp = state.pressure - state.w_p;
                dV = V - state.w;
                dVdp = bsxfun(@rdivide, dV, dp);
                dVdp(~isfinite(dVdp)) = 0;
            else
                dVdp = 0*V; 
            end
        case 'eos'
            % Already computed if valid
        case 'none'
            dVdp = zeros(size(V));
        case 'analytical'
            % Already computed if valid
        otherwise
            error('Unknown derivative strategy %s', dstrategy);
    end
end

function [V, dVdp] = getPartialVolumesInternal(model, state, opt, subs, computeDerivatives)
    pv = model.operators.pv(subs);
    V = getTwoPhaseVolumes(model, state, pv);
    if computeDerivatives && nargout > 1
        state_perturb = state;
        dp = opt.pressureEpsilon;
        state_perturb.pressure = state_perturb.pressure - dp;
        state_perturb = model.computeFlash(state_perturb, inf);
        
        V_lower = getTwoPhaseVolumes(model, state_perturb, pv);
        dVdp = (V - V_lower)./dp;
    else
        dVdp = 0*V;
    end

    switch lower(opt.units)
        case 'mass'
            V = bsxfun(@rdivide, V, model.EOSModel.fluid.molarMass);
        case 'moles'
            % Nothing to do here
        otherwise
            error('Unknown system %s', opt.type);
    end
end

function state = makeSubstate(state, flag)
    flds = getCompCellFields();
    for i = 1:numel(flds)
        f = flds{i};
        if isfield(state, f)
            state.(f) = state.(f)(flag, :);
        end
    end
end

function V = getTwoPhaseVolumes(model, state, pv)
    [x, y, p, T, Z_L, Z_V, L] = ...
        model.getProps(state, ...
        'x', 'y', 'pressure', 'T', 'Z_L', 'Z_V', 'L');
    z = state.components;
    ncomp = size(x, 2);
    if isempty(pv)
        V = zeros(0, ncomp);
        return
    end
    
    V = 1 - L;
    
    sL = state.s(:, 1+model.water);
    sV = state.s(:, 2+model.water);
    
    rhoL = model.EOSModel.PropertyModel.computeMolarDensity(p, x, Z_L, T, true);
    rhoV = model.EOSModel.PropertyModel.computeMolarDensity(p, y, Z_V, T, false);

    N_L = pv.*rhoL.*sL;
    N_V = pv.*rhoV.*sV;
    N_T = N_L.*L + N_V.*V;
    
    z = expandMatrixToCell(z);
    N = z;
    for i = 1:ncomp
        N{i} = z{i}.*N_T;
    end
    
    [N{:}] = initVariablesADI(N{:});
    
    Nt = 0;
    for i = 1:ncomp
        Nt = Nt + N{i};
    end
    
    z = N;
    for i = 1:ncomp
        z{i} = z{i}./Nt;
    end
    
    state.eos.packed = model.EOSModel.getPropertiesFastAD(state.pressure, state.T, state.x, state.y);
    [x, y, L, packed] = model.EOSModel.getPhaseFractionAsADI(state, p, T, z);
    
    [Z_L, Z_V] = model.EOSModel.getCompressibility(state, p, T, x, y, z);
    
    rhoL = model.EOSModel.PropertyModel.computeMolarDensity(p, x, Z_L, T, true);
    rhoV = model.EOSModel.PropertyModel.computeMolarDensity(p, y, Z_V, T, false);
    Vt = Nt.*(L./rhoL + (1-L)./rhoV);
    
    V = zeros(numel(pv), ncomp);
    for i = 1:ncomp
        V(:, i) = diag(Vt.jac{i});
    end
end

function [V, dVdp] = getSinglePhaseVolumes(model, state, liquid, vapor, computeDerivatives)
    singlePhase = liquid | vapor;
    if ~any(singlePhase)
        [V, dVdp] = deal(zeros(0, size(state.components, 2)));
        return
    end
    z = state.components(singlePhase, :);
    p = state.pressure(singlePhase);
    Z_L = state.Z_L(singlePhase);
    Z_V = state.Z_V(singlePhase);
    T = state.T(singlePhase);
    if computeDerivatives
        p = initVariablesAD_diagonal(p);
        Z_L = double2GenericAD(Z_L, p);
        Z_V = double2GenericAD(Z_V, p);
        eos = model.EOSModel;
        zc = expandMatrixToCell(z);
        [Si_L, Si_V, A_L, A_V, B_L, B_V] =...
            eos.getMixtureFugacityCoefficients(p, T, zc, zc, eos.fluid.acentricFactors);
        Z_L = eos.setZDerivatives(Z_L, A_L, B_L);
        Z_V = eos.setZDerivatives(Z_V, A_V, B_V);
    end
    
    liquid = liquid(singlePhase);
    vapor = vapor(singlePhase);

    Z = Z_L.*liquid + Z_V.*vapor;
    rho = model.EOSModel.PropertyModel.computeMolarDensity(p, z, Z, T, nan);
    
    V = 1./rho;
    ncomp = model.EOSModel.fluid.getNumberOfComponents();
    if computeDerivatives
        dVdp = V.jac{1}.diagonal;
        dVdp = repmat(dVdp, 1, ncomp);
        V = value(V);
    end
    V = repmat(V, 1, ncomp);
end

function flds = getCompCellFields()
    flds = {'pressure', 's', 'x', 'y', 'components', 'K', ...
            'mob', 'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V', 'w_p', 'w'};
end

function w = getWeights(problem, subs)
    acc = problem.accumulationTerms;
    state = problem.state;
    [ncell, ncomp] = size(problem.state.components);
    hasWater = size(state.s, 2) == 3;
    if hasWater
        ncomp = ncomp + 1;
    end
    c = combineEquations(acc{:});
    if isnumeric(c)
        if nargin > 1
            ncell = sum(subs);
        end
        w = ones(ncell, ncomp);
        return;
    end
    J = c.jac{1};

    if ~isempty(problem.reorder)
        J = J(:, problem.reorder);
    end

    J(:, problem.wellVarIndices) = [];

    ndof = ncell*ncomp;
    [B, C, D, E] = getBlocks(J, ndof);
    [L, U] = lu(E);
    A = B - C*(U\(L\D));
%             A = B - C*(E\D);
    b = zeros(ndof, 1);
    b(1:ncell) = 1/barsa;

    w = (A')\b;
    w = reshape(w, [], ncomp);
    w = bsxfun(@rdivide, w, sum(abs(w), 2));
    w = bsxfun(@rdivide, w, sum(state.rho.*state.s, 2));
    if nargin > 1
        w = w(subs, :);
    end
end

function [B, C, D, E] = getBlocks(J, ndof)
    start = 1:ndof;

    if 0
        stop = (ndof+1):size(J, 2);
        B = J(start, start);
        C = J(start, stop);
        D = J(stop, start);
        E = J(stop, stop);
    else
       [ix, jx, vx] = find(J);
       n = size(J, 2);
       keep = false(n, 1);
       keep(start) = true;
       nk = ndof;

       keepRow = keep(ix);
       keepCol = keep(jx);
       kb = keepRow & keepCol;
       B = sparse(ix(kb), jx(kb), vx(kb), nk, nk);

       kc = keepRow & ~keepCol;
       C = sparse(ix(kc), jx(kc) - nk, vx(kc), nk, n - nk);

       kd = ~keepRow & keepCol;
       D = sparse(ix(kd) - nk, jx(kd), vx(kd), n - nk, nk);

       ke = ~keepRow & ~keepCol;
       E = sparse(ix(ke) - nk, jx(ke) - nk, vx(ke), n - nk, n - nk);
    end

end
