classdef NatVarsShaleModel < NaturalVariablesCompositionalModel
    % Natural variables model for compositional problems
    %
    % SYNOPSIS:
    %   model = NaturalVariablesCompositionalModel(G, rock, fluid, compFluid)
    %
    % DESCRIPTION:
    %   The natural variables model relies on separate primary variables
    %   for saturations in the multiphase region. This makes it possible to
    %   perform certain heuristics to ensure convergence (e.g. saturation
    %   chopping) since there is a weaker correspondence between
    %   saturations and the compositions. However, this formulation
    %   traverses the phase boundary as a number of discrete steps and can
    %   take more nonlinear iterations for some problems.
    %
    % PARAMETERS:
    %   G         - Grid structure
    %   rock      - Rock structure for the reservoir
    %   fluid     - The flow fluid, containing relative permeabilities,
    %               surface densities and flow properties for the
    %               aqueous/water phase (if present)
    %   compFluid - CompositionalMixture instance describing the species
    %               present.
    %
    % RETURNS:
    %  model - Initialized class instance
    %
    % SEE ALSO:
    %   `NaturalVariablesCompositionalModel, ThreePhaseCompositionalModel`, `OverallCompositionCompositionalModel`
    
    methods
        function model = NatVarsShaleModel(G, rock, fluid, compFluid, varargin)
            model = model@NaturalVariablesCompositionalModel(G, rock, fluid, compFluid, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = eqnsShaleNaturalVars(state0, state, model, dt, ...
                            drivingForces, varargin{:});
        end
    end
end

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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
