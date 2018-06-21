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
    
    getRegCell = @(x) repmat(opt.cells(1), size(double(x)));
    
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
        
        rho{ix} = @(p, z) getOilDensity(model, p, z, opt.rs, 'cellInx', getRegCell(p));
        PC{ix} = @(S) 0*S;
    end
    
    if model.gas
        ix = model.getPhaseIndex('G');
        rho{ix} = @(p, z) getGasDensity(model, p, z, opt.rv, 'cellInx', getRegCell(p));
        pc_sign(ix) = 1;
        if isfield(model.fluid, 'pcOG')
            PC{ix} = @(S) model.fluid.pcOG(S, 'cellInx', getRegCell(S));
        else
            PC{ix} = @(S) 0*S;
        end
    end
    if model.oil
        ref_index = model.getPhaseIndex('O');
    else
        ref_index = 1;
    end
    
    [s_min, s_max] = getMinMaxPhaseSaturations(model, opt.cells);
    
    region = getInitializationRegionsBase(model, rho, contacts, ...
        'rho',              rho, ...
        'cells',            opt.cells, ...
        'reference_index',  ref_index, ...
        'pc_sign',          pc_sign, ...
        's_min',            s_min, ...
        's_max',            s_max, ...
        'pc_functions',     PC, ...
        args{:});
    region.rs = opt.rs;
    region.rv = opt.rv;
end

function rhoO = getOilDensity(model, p, z, rs, varargin)
    f = model.fluid;
    if model.disgas
        assert(~isempty(rs), 'Disgas model must have rs specified as optional argument');
        if isa(rs, 'function_handle')
            rs = rs(p, z);
        end
        rsSat = f.rsSat(p, varargin{:});
        rs = min(rs, rsSat);
        rhoO = f.bO(p, rs, rs >= rsSat, varargin{:}).*(f.rhoOS + rs.*f.rhoGS);
    else
        rhoO = f.bO(p, varargin{:}).*f.rhoOS;
    end
end

function rhoG = getGasDensity(model, p, z, rv, varargin)
    f = model.fluid;
    if model.vapoil
        assert(~isempty(rv), 'Vapoil model must have rv specified as optional argument');
        if isa(rv, 'function_handle')
            rv = rv(p, z);
        end
        rvSat = f.rvSat(p, varargin{:});
        rv = min(rv, rvSat);
        rhoG = f.bG(p, rv, rv >= rvSat, varargin{:}).*(rv.*f.rhoOS + f.rhoGS);
    else
        rhoG = f.bG(p, varargin{:}).*f.rhoGS;
    end
end