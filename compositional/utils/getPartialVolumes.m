function [V, dVdp] = getPartialVolumes(model, state, varargin)
    opt = struct('type', 'mass', ...
                 'simple_singlephase', true);
    opt = merge_options(opt, varargin{:});
    % Get partial molar volumes
    V = getPartialVolumesInternal(model, state, opt);
    if nargout > 1
        state1 = state;
        dp = mean(state.pressure)*1e-6;
        state1.pressure = state1.pressure - dp;
        state1 = model.computeFlash(state1, inf);
        
        V_lower = getPartialVolumesInternal(model, state1, opt);
        dVdp = (V - V_lower)./dp;
    end
end
function V = getPartialVolumesInternal(model, state, opt)
    nc = model.G.cells.num;
    if opt.simple_singlephase
        [pureLiquid, pureVapor, twoPhase] = model.getFlag(state);
    else
        pureLiquid = false(nc, 1);
        pureVapor = false(nc, 1);
        twoPhase = true(nc, 1);
    end
    ncomp = model.EOSModel.fluid.getNumberOfComponents();
    V = zeros(nc, ncomp);
    V(~twoPhase, :) = getSinglePhaseVolumes(model, state, pureLiquid, pureVapor);
    V(twoPhase, :) = getTwoPhaseVolumes(model, state, twoPhase);
    switch lower(opt.type)
        case 'mass'
            V = bsxfun(@rdivide, V, model.EOSModel.fluid.molarMass);
        case 'moles'
        otherwise
            error('Unknown system %s', opt.type);
    end
end

function V = getTwoPhaseVolumes(model, state, twoPhase)
    flds = getCompCellFields();
    for i = 1:numel(flds)
        f = flds{i};
        if isfield(state, f)
            state.(f) = state.(f)(twoPhase, :);
        end
    end

    pv = model.operators.pv(twoPhase);
    
    [x, y, p, T, Z_L, Z_V, L] = ...
        model.getProps(state, ...
        'x', 'y', 'pressure', 'T', 'Z_L', 'Z_V', 'L');
    z = state.components;
    ncomp = size(x, 2);
    
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
%     
%     V = 1./rho;
    Vt = Nt.*(L./rhoL + (1-L)./rhoV);
    
    V = zeros(nnz(twoPhase), ncomp);
    for i = 1:ncomp
        V(:, i) = diag(Vt.jac{i});
    end
end

function V = getSinglePhaseVolumes(model, state, liquid, vapor)
    singlePhase = liquid | vapor;
    z = state.components(singlePhase, :);
    p = state.pressure(singlePhase);
    T = state.T(singlePhase);
    
    liquid = liquid(singlePhase);
    vapor = vapor(singlePhase);
    
    
    Z = state.Z_L(singlePhase).*liquid + state.Z_V(singlePhase).*vapor;
    rho = model.EOSModel.PropertyModel.computeMolarDensity(p, z, Z, T, nan);
    
    V = 1./rho;
    ncomp = model.EOSModel.fluid.getNumberOfComponents();
    V = repmat(V, 1, ncomp);
end

function flds = getCompCellFields()
    flds = {'pressure', 's', 'x', 'y', 'components', 'K', ...
            'mob', 'rho', 'flag', 'L', 'T', 'Z_L', 'Z_V'};
end

