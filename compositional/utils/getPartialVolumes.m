function [V, dVdp] = getPartialVolumes(model, state, varargin)
    opt = struct('type', 'mass', ...
                 'pressure_perturb', mean(state.pressure)*1e-3, ...
                 'simple_singlephase', false);
    opt = merge_options(opt, varargin{:});
    computeDerivatives = nargout > 1;
    
    % Get partial molar volumes
    [V, dVdp] = getPartialVolumesInternal(model, state, opt, computeDerivatives);
end

function [V, dVdp] = getPartialVolumesInternal(model, state, opt, computeDerivatives)
    nc = model.G.cells.num;
    if opt.simple_singlephase
        [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
    else
        pureLiquid = false(nc, 1);
        pureVapor = false(nc, 1);
        twoPhase = true(nc, 1);
    end
    ncomp = model.EOSModel.fluid.getNumberOfComponents();
    [V, dVdp] = deal(zeros(nc, ncomp));
    % Single phase regime can get derivatives from eos in regular wau
    [V(~twoPhase, :), dVdp(~twoPhase, :)] = getSinglePhaseVolumes(model, state, pureLiquid, pureVapor, computeDerivatives);
    
    % Numerical perturbation w/flash to get two-phase region weights
    substate = makeSubstate(state, twoPhase);
    pv = model.operators.pv(twoPhase);
    V(twoPhase, :) = getTwoPhaseVolumes(model, substate, pv);
    if computeDerivatives
        state_perturb = substate;
        dp = opt.pressure_perturb;
        state_perturb.pressure = state_perturb.pressure - dp;
        state_perturb = model.computeFlash(state_perturb, inf);
        
        V_lower = getTwoPhaseVolumes(model, state_perturb, pv);
        dVdp(twoPhase, :) = (V(twoPhase, :) - V_lower)./dp;
    end
    
    switch lower(opt.type)
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
%     state = makeSubstate(state, twoPhase);
%     pv = model.operators.pv(twoPhase);
    
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
        Z_L = double2NewAD(Z_L, p);
        Z_V = double2NewAD(Z_V, p);
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
        V = double(V);
    end
    V = repmat(V, 1, ncomp);
end

function flds = getCompCellFields()
    flds = {'pressure', 's', 'x', 'y', 'components', 'K', ...
            'mob', 'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V'};
end

