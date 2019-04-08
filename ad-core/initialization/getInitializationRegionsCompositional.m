function region = getInitializationRegionsCompositional(model, contacts, varargin)
    opt = struct('cells',   (1:model.G.cells.num)', ...
                 'T',       303.15, ...
                 'x',       [], ...
                 'y',       []);
    [opt, args] = merge_options(opt, varargin{:});
    x = opt.x;
    y = opt.y;
    T = opt.T;
    actPh = model.getActivePhases();
    nPh = sum(actPh);

    rho = cell(1, nPh);
    PC = cell(1, nPh);
    pc_sign = ones(1, nPh);
    
    
    
    f = model.fluid;
    [satnum, pvtnum] = deal(1);
    if isfield(model.rock, 'regions')
        if isfield(model.rock.regions, 'saturation')
            satnum = model.rock.regions.saturation(cc);
        end
        if isfield(model.rock.regions, 'pvt')
            pvtnum = model.rock.regions.pvt(cc);
        end
    end
    if model.water
        ix = model.getPhaseIndex('W');
        
        bW = getFunction(f, 'bW', pvtnum);
        rho{ix} = @(p, z) bW(p).*f.rhoWS(pvtnum);
        pc_sign(ix) = -1;
        if isfield(model.fluid, 'pcOW')
            pcOW = getFunction(f, 'pcOW', satnum);
            PC{ix} = @(S) pcOW(S);
        else
            PC{ix} = @(S) 0*S;
        end
    end
    
    if model.oil
        ix = model.getPhaseIndex('O');
        rho{ix} = @(p, z) getDensity(model, p, T, z, x, true);
        PC{ix} = @(S) 0*S;
    end
    
    if model.gas
        ix = model.getPhaseIndex('G');
        rho{ix} = @(p, z) getDensity(model, p, T, z, y, false);
        pc_sign(ix) = 1;
        if isfield(model.fluid, 'pcOG')
            pcOG = getFunction(f, 'pcOG', satnum);
            PC{ix} = @(S) pcOG(S);
        elseif ~model.oil && isfield(model.fluid, 'pcWG')
            pcWG = getFunction(f, 'pcWG', satnum);
            PC{ix} = @(S) pcWG(S);
        else
            PC{ix} = @(S) 0*S;
        end
    end
    ref_index = model.getPhaseIndex('O');
    
    [s_min, s_max] = getMinMaxPhaseSaturations(model, satnum, opt.cells);
    
    region = getInitializationRegionsBase(model, rho, contacts, ...
        'rho',              rho, ...
        'cells',            opt.cells, ...
        'reference_index',  ref_index, ...
        'pc_sign',          pc_sign, ...
        's_min',            s_min, ...
        's_max',            s_max, ...
        'pc_functions',     PC, ...
        args{:});
end

function rho = getDensity(model, p, T, z, x, isLiquid, varargin)
    eos = model.EOSModel;
    if isa(model.EOSModel, 'EquilibriumConstantModel')
        Z = nan;
    else
        [A_ij, Bi] = eos.getMixingParameters(p, T, eos.fluid.acentricFactors, iscell(x));
        [Si, A, B] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
        Z = model.EOSModel.computeCompressibilityZ(p, x, A, B, Si, Bi, isLiquid);
    end
    rho = model.EOSModel.PropertyModel.computeDensity(p, x, Z, T, isLiquid);
end

function f = getFunction(fluid, fld, reg)
    f = fluid.(fld);
    if iscell(f)
        f = f{reg};
    end 
end