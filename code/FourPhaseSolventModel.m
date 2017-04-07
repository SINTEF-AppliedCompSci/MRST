classdef FourPhaseSolventModel < ReservoirModel
% Three phase with optional dissolved gas and vaporized oil
properties
   solvent
end

methods
    function model = FourPhaseSolventModel(G, rock, fluid, varargin)
        model = model@ReservoirModel(G, rock, fluid, varargin{:});

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        % All phases are present
        model.water = true;
        model.oil = true;
        model.gas = true;
        model.solvent = true;
        model.saturationVarNames = {'sw', 'so', 'sg', 'ss'};

        model = merge_options(model, varargin{:});

    end
    
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        switch(lower(name))
            case {'solvent', 'ss'}
                index = 4;
                fn = 's';
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ReservoirModel(model, name);
        end
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsFourPhaseSolvent(state0, state, ...
                model, dt, drivingForces, varargin{:});

    end

    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@ReservoirModel(model, state);
    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces)
        saturations = lower(model.saturationVarNames);
        wi = strcmpi(saturations, 'sw');
        oi = strcmpi(saturations, 'so');
        gi = strcmpi(saturations, 'sg');

        % Parent class handles almost everything for us
        [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);

        % Handle the directly assigned values (i.e. can be deduced directly from
        % the well controls. This is black oil specific.
        W = drivingForces.W;
        state.wellSol = assignWellValuesFromControl(model, state.wellSol, W, wi, oi, gi);
    end
    
     function [isActive, phInd] = getActivePhases(model)
        % Get active flag for the canonical phase ordering (water, oil
        % gas as on/off flags).
        isActive = [model.water, model.oil, model.gas, model.solvent];
        if nargout > 1
            phInd = find(isActive);
        end
    end
    
end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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