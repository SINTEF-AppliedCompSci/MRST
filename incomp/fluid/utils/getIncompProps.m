function varargout = getIncompProps(state, fluid)
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

    varargout = cell(1, nargout);
    if isfield(fluid, 'properties')
        [varargout{:}] = getPropsLegacy(state, fluid);
    else
        [varargout{:}] = getPropsAD(state, fluid);
    end
end

function [rho, kr, mu, dkr] = getPropsLegacy(state, fluid)
    [mu, rho] = fluid.properties(state);
    s         = fluid.saturation(state);
    if nargout < 4
        kr = fluid.relperm(s, state);
    else
        [kr, dkr] = fluid.relperm(s, state);
    end
end

function [rho, kr, mu, dkr] = getPropsAD(state, fluid)
    if isfield(state, 'pressure')
        p = state.pressure;
        p_avg = mean(p);
    else
        p_avg = 1*atm;
    end
    
    rho = [fluid.rhoWS*fluid.bW(p_avg), ...
           fluid.rhoOS*fluid.bO(p_avg)];
    mu = [fluid.muW(p), ...
           fluid.muO(p)];
    getDer = nargout > 3;
    s = state.s(:, 1);
    if getDer
        s = initVariablesAD_diagonal(s);
    end
    
    krW = fluid.krW(s);
    so = 1 - s;
    if isfield(fluid, 'krOW')
        krO = fluid.krOW(so);
    else
        krO = fluid.krO(so);
    end
    kr = [value(krW), value(krO)];
    if getDer
        dkr = [krW.jac{1}.diagonal, krO.jac{1}.diagonal]; 
    end
end
