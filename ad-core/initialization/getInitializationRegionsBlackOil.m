function region = getInitializationRegionsBlackOil(model, contacts, varargin)
    opt = struct('cells', (1:model.G.cells.num)',...
                 'rs', [], ...
                 'rv', []);
    [opt, args] = merge_options(opt, varargin{:});
    
    actPh = model.getActivePhases();
    nPh = sum(actPh);
    
    rho = cell(1, nPh);
    PC = cell(1, nPh);
    pc_sign = ones(1, nPh);
    
    cc = opt.cells(1);
    
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
        
        rho{ix} = @(p, z) getOilDensity(model, p, z, opt.rs, pvtnum);
        PC{ix} = @(S) 0*S;
    end
    
    if model.gas
        ix = model.getPhaseIndex('G');
        rho{ix} = @(p, z) getGasDensity(model, p, z, opt.rv, pvtnum);
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
    if model.oil
        ref_index = model.getPhaseIndex('O');
    else
        ref_index = 1;
    end
    
    
    [s_min, s_max] = getMinMaxPhaseSaturations(model, satnum);
    
    region = getInitializationRegionsBase(model, rho, contacts, ...
        'rho',              rho, ...
        'cells',            opt.cells, ...
        'reference_index',  ref_index, ...
        'pc_sign',          pc_sign, ...
        's_min',            s_min, ...
        's_max',            s_max, ...
        'saturation_region',           satnum, ...
        'pvt_region',              pvtnum, ...
        'pc_functions',     PC, ...
        args{:});
    region.rs = opt.rs;
    region.rv = opt.rv;
end

function rhoO = getOilDensity(model, p, z, rs, pvtnum)
    f = model.fluid;
    bO = getFunction(f, 'bO', pvtnum);
    if model.disgas
        assert(~isempty(rs), 'Disgas model must have rs specified as optional argument');
        if isa(rs, 'function_handle')
            rs = rs(p, z);
        end
        rsSatF = getFunction(f, 'rsSat', pvtnum);
        rsSat = rsSatF(p);
        rs = min(rs, rsSat);
        rhoO = bO(p, min(rs, rsSat), rs >= rsSat).*(f.rhoOS(pvtnum) + rs.*f.rhoGS(pvtnum));
    else
        rhoO = bO(p).*f.rhoOS(pvtnum);
    end
end

function rhoG = getGasDensity(model, p, z, rv, pvtnum)
    f = model.fluid;
    bG = getFunction(f, 'bG', pvtnum);
    if model.vapoil
        assert(~isempty(rv), 'Vapoil model must have rv specified as optional argument');
        if isa(rv, 'function_handle')
            rv = rv(p, z);
        end
        rvSatF = getFunction(f, 'rvSat', pvtnum);
        rvSat = rvSatF(p);
        rv = min(rv, rvSat);
        rhoG = bG(p, min(rv, rvSat), rv >= rvSat).*(rv.*f.rhoOS(pvtnum) + f.rhoGS(pvtnum));
    else
        rhoG = bG(p).*f.rhoGS(pvtnum);
    end
end

function f = getFunction(fluid, fld, reg)
    f = fluid.(fld);
    if iscell(f)
        f = f{reg};
    end 
end