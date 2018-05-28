function region = getInitializationRegionsCompositional(model, cells, datum_p, datum_z, contacts, contacts_pc, x, y, T)
    actPh = model.getActivePhases();
    nPh = sum(actPh);
    
    if ischar(cells)
        getRegCell = @(x) repmat(1, size(double(x)));
    else
        getRegCell = @(x) repmat(cells(1), size(double(x)));
    end
    
    rho = cell(1, nPh);
    PC = cell(1, nPh);
    pc_sign = ones(1, nPh);
    
    
    
    f = model.fluid;
    if model.water
        ix = model.getPhaseIndex('W');
        
        rho{ix} = @(p, z) f.bW(p, 'cellInx', getRegCell(p)).*f.rhoWS;
        pc_sign(ix) = -1;
        if isfield(model.fluid, 'pcOW')
            PC{ix} = @(S) model.fluid.pcOW(S, 'cellInx', getRegCell(S));
        else
            PC{ix} = @(S) 0*S;
        end
    end
    
    if model.oil
        ix = model.getPhaseIndex('O');
        rho{ix} = @(p, z) getDensity(model, p, T, z, x, true, 'cellInx', getRegCell(p));
        PC{ix} = @(S) 0*S;
    end
    
    if model.gas
        ix = model.getPhaseIndex('G');
        rho{ix} = @(p, z) getDensity(model, p, T, z, y, false, 'cellInx', getRegCell(p));
        pc_sign(ix) = 1;
        if isfield(model.fluid, 'pcOG')
            PC{ix} = @(S) model.fluid.pcOG(S, 'cellInx', getRegCell(S));
        else
            PC{ix} = @(S) 0*S;
        end
    end
    ref_index = model.getPhaseIndex('O');
    
    [s_min, s_max] = getMinMaxPhaseSaturations(model, cells);
    
    region = getInitializationRegionsBase(model, cells, ref_index, rho, datum_p, datum_z, contacts, ...
        'contacts_pc', contacts_pc, ...
        'pc_sign',      pc_sign, ...
        's_min',        s_min, ....
        's_max',        s_max, ....
        'pc_functions', PC);
end

function rho = getDensity(model, p, T, z, x, isLiquid, varargin)
    eos = model.EOSModel;
    if isa(model.EOSModel, 'EquilibriumConstantModel')
        Z = nan;
    else
        [A_ij, Bi] = eos.getMixingParameters(p, T, eos.fluid.acentricFactors, iscell(x));
        [Si, A, B] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
        Z = model.EOSModel.computeCompressibilityZ(p, x, A, B, Si, Bi);
    end
    rho = model.EOSModel.PropertyModel.computeDensity(p, x, Z, T, isLiquid);
end