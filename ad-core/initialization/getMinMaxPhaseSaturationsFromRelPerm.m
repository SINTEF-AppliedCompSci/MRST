function [s_min, s_max] = getMinMaxPhaseSaturationsFromRelPerm(model, tol, reg)
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

    if nargin < 2
        tol = 1e-6;
    end
    relPerm = model.FlowPropertyFunctions.RelativePermeability;
    getKr = @(name, s) relPerm.evaluateFunctionSingleRegionWithArguments(model.fluid.(['kr', name]), reg, s);
    nph = sum(model.getActivePhases());
    s = (0:tol:1)';
    
    if nargin < 3
        reg = 1;
    end
    
    s_min = zeros(1, nph);
    s_max = ones(1, nph);
    
    if model.water
        krw = getKr('W', s);
        
        ix = model.getPhaseIndex('W');
        s_min(ix) = getMinSat(s, krw);
        s_max(ix) = getMaxSat(s, krw);
    end
    
    if model.oil
        if isfield(model.fluid, 'krO')
            kro = getKr('O', s);
        else
            if model.water
                krow = getKr('OW', s);
            else
                krow = 1;
            end
            
            if model.gas
                if isfield(model.fluid, 'krO')
                    krog = getKr('O', s);
                else
                    krog = getKr('OG', s);
                end
            else
                krog = 1;
            end
            kro = min(krow, krog);
        end
        
        ix = model.getPhaseIndex('O');
        s_min(ix) = getMinSat(s, kro);
        s_max(ix) = getMaxSat(s, kro);
    end
    
    if model.gas
        krg = getKr('G', s);
        
        ix = model.getPhaseIndex('G');        
        s_min(ix) = getMinSat(s, krg);
        s_max(ix) = getMaxSat(s, krg);
    end
end

function s_min = getMinSat(s, kr)
    sub = find(kr == 0, 1, 'last');
    s_min = s(sub);
end

function s_max = getMaxSat(s, kr)
    sub = find(kr == kr(end), 1, 'first');
    s_max = s(sub);
end
