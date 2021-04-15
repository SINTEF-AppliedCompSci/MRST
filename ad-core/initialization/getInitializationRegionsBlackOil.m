function region = getInitializationRegionsBlackOil(model, contacts, varargin)
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

    opt = struct('cells', (1:model.G.cells.num)',...
                 'rs', [], ...
                 'rv', []);
    [opt, args] = merge_options(opt, varargin{:});
    % Ensure that groups are set up before building init regions
    model = model.setupStateFunctionGroupings();

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
    [pc, pc_sign, pc_scale] = getEquilPC(model, satnum, opt.cells);
    if model.water
        ix = model.getPhaseIndex('W');
        bW = getFunction(f, 'bW', pvtnum);
        rho{ix} = @(p, z) bW(p).*f.rhoWS(pvtnum);
    end
    
    if model.oil
        ix = model.getPhaseIndex('O');
        rho{ix} = @(p, z) getOilDensity(model, p, z, opt.rs, pvtnum);
    end
    
    if model.gas
        ix = model.getPhaseIndex('G');
        rho{ix} = @(p, z) getGasDensity(model, p, z, opt.rv, pvtnum);
    end
    if model.oil
        ref_index = model.getPhaseIndex('O');
    else
        ref_index = 1;
    end
    [s_min, s_max] = getMinMaxPhaseSaturations(model, satnum, opt.cells);
    region = getInitializationRegionsBase(model, rho, contacts, ...
        'rho',               rho, ...
        'cells',             opt.cells, ...
        'reference_index',   ref_index, ...
        'pc_sign',           pc_sign, ...
        'pc_scale',          pc_scale, ...
        's_min',             s_min, ...
        's_max',             s_max, ...
        'saturation_region', satnum, ...
        'pvt_region',        pvtnum, ...
        'pc_functions',      pc, ...
        args{:});
    region.rs = opt.rs;
    region.rv = opt.rv;
end

function rhoO = getOilDensity(model, p, z, rs, pvtnum)
    f = model.fluid;
    bO = getFunction(f, 'bO', pvtnum);
    if model.disgas
        p = max(p, model.minimumPressure);
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
        p = max(p, model.minimumPressure);
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
