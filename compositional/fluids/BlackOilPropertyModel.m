classdef BlackOilPropertyModel < PropertyModel
    % Black-oil property model for compositional models
    properties
        viscosityFunctions
        densityFunctions
        checkStabilityFunction
    end
    
    methods
        function model = BlackOilPropertyModel(rho, mu, satfn, varargin)
            model = model@PropertyModel(varargin{:});
            model.densityFunctions = rho;
            model.viscosityFunctions = mu;
            model.checkStabilityFunction = satfn;
        end
        function rho = computeDensity(model, eos, p, x, Z, T, isLiquid)
            if isLiquid
                ix = 1;
            else
                ix = 2;
            end
            x = expandMatrixToCell(x);
            rho = model.densityFunctions{ix}(p, T, x);
        end
        
        function rho = computeMolarDensity(model, eos, p, x, Z, T, isLiquid)
            x = expandMatrixToCell(x);
            rho = model.computeDensity(eos, p, x, Z, T, isLiquid);
        end
        
        function mu = computeViscosity(model, eos, P, x, Z, T, isLiquid)
            if isLiquid
                ix = 1;
            else
                ix = 2;
            end
            x = expandMatrixToCell(x);
            mu = model.viscosityFunctions{ix}(P, T, x);
        end
    end
end

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

