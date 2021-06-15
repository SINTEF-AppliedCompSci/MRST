function region = getInitializationRegionsCompositional(model, contacts, varargin)
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

    opt = struct('cells',   (1:model.G.cells.num)', ...
                 'T',       [], ...
                 'x',       [], ...
                 'y',       []);
    [opt, args] = merge_options(opt, varargin{:});
    x = makeFunction(opt.x);
    y = makeFunction(opt.y);
    if isempty(opt.T)
        opt.T = 303.15;
    end
    T = makeFunction(opt.T);
    actPh = model.getActivePhases();
    nPh = sum(actPh);

    rho = cell(1, nPh);

    f = model.fluid;
    [satnum, pvtnum] = deal(1);
    if isfield(model.rock, 'regions')
        if isfield(model.rock.regions, 'pvt')
            pvtnum = model.rock.regions.pvt(opt.cells);
        end
        if isfield(model.rock.regions, 'saturation')
            satnum = model.rock.regions.saturation(opt.cells);
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
        rho{ix} = @(p, z) getDensity(model, p, T(p, z), z, x(p, z), true);
    end
    if model.gas
        ix = model.getPhaseIndex('G');
        rho{ix} = @(p, z) getDensity(model, p, T(p, z), z, y(p, z), false);
    end
    ref_index = model.getPhaseIndex('O');
    
    [s_min, s_max] = getMinMaxPhaseSaturations(model, satnum, opt.cells);
    
    region = getInitializationRegionsBase(model, rho, contacts, ...
        'rho',              rho, ...
        'cells',            opt.cells, ...
        'reference_index',  ref_index, ...
        'pc_sign',          pc_sign, ...
        'pc_scale',         pc_scale, ...
        's_min',            s_min, ...
        's_max',            s_max, ...
        'pc_functions',     pc, ...
        args{:});
end

function rho = getDensity(model, p, T, z, x, isLiquid, varargin)
    eos = model.EOSModel;
    if isa(model.EOSModel, 'EquilibriumConstantModel')
        Z = nan;
    else
        [A_ij, Bi] = eos.getMixingParameters(p, T, eos.CompositionalMixture.acentricFactors, iscell(x));
        [Si, A, B] = eos.getPhaseMixCoefficients(x, A_ij, Bi);
        Z = eos.computeCompressibilityZ(p, x, A, B, Si, Bi, isLiquid);
    end
    rho = eos.PropertyModel.computeDensity(eos, p, x, Z, T, isLiquid);
end

function f = makeFunction(f)
    if ~isa(f, 'function_handle')
        f = @(p, z) 0*z + f;
    end
end

function f = getFunction(fluid, fld, reg)
    f = fluid.(fld);
    if iscell(f)
        f = f{reg};
    end 
end
