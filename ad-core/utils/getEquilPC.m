function [PC, pc_sign, pc_scalers] = getEquilPC(model, satnum, cells)
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

    actPh = model.getActivePhases();
    nPh = sum(actPh);
    f = model.fluid;

    nc = numel(cells);
    PC = cell(1, nPh);
    pc_sign = ones(1, nPh);
    pc_scalers = ones(nc, nPh);

    pcimpl = model.FlowPropertyFunctions.CapillaryPressure;
    
    if model.water
        ix = model.getPhaseIndex('W');
        
        pc_sign(ix) = -1;
        if isfield(model.fluid, 'pcOW')
            pcOW = getFunction(f, 'pcOW', satnum);
            PC{ix} = @(S) pcOW(S);
            if pcimpl.hasJFunctionScaler('OW')
                ratio = pcimpl.getJFunctionStaticRatio(model);
                pc_scalers(:, ix) = ratio(cells).*pcimpl.getSurfaceTension('OW');
            end
        else
            PC{ix} = @(S) 0*S;
        end
    end
    if model.oil
        ix = model.getPhaseIndex('O');
        PC{ix} = @(S) 0*S;
    end
    if model.gas
        ix = model.getPhaseIndex('G');
        pc_sign(ix) = 1;
        if isfield(model.fluid, 'pcOG')
            pcOG = getFunction(f, 'pcOG', satnum);
            PC{ix} = @(S) pcOG(S);
            if pcimpl.hasJFunctionScaler('OG')
                ratio = pcimpl.getJFunctionStaticRatio(model);
                pc_scalers(:, ix) = ratio(cells).*pcimpl.getSurfaceTension('OG');
            end
        elseif ~model.oil && isfield(model.fluid, 'pcWG')
            pcWG = getFunction(f, 'pcWG', satnum);
            PC{ix} = @(S) pcWG(S);
        else
            PC{ix} = @(S) 0*S;
        end
    end
end

function f = getFunction(fluid, fld, reg)
    f = fluid.(fld);
    if iscell(f)
        f = f{reg};
    end 
end
