classdef DualContWaterFluidModel < DualContinuumWaterModel
%
% SYNOPSIS:
%   model = DualContWaterFluidModel(G, rock, fluid, varargin)
%
% DESCRIPTION: 
%   This model is for the fluid part of a fully coupled dual continuum 
%   poroelastic simulation. It is derived from the fluid model DualContWaterFluidModel 
%   and a few functionalities that are needed for the coupled solver are added.
%
% PARAMETERS:
%   G     - Simulation grid.
%   rock  - Rock cell for fracture / matrix
%   fluid - Fluid cell for fracture / matrix
%
% RETURNS:
%   class instance
%
% EXAMPLE: 
%
% SEE ALSO: DualContMechFluidModel 
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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

    properties
        primaryVarNames;
    end

    methods
        function model = DualContWaterFluidModel(G, rock, fluid, varargin)
            model = model@DualContinuumWaterModel(G, rock, fluid);
            model.primaryVarNames = {'pressure', 'pressure_matrix'}; % well variables
                                                                     % not included
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            error(['The DualContWaterFluidModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from DualContMechFluidModel.'])
        end

        function varnames = getAllPrimaryVariables(model, state)
        % List all the primary variables that are recognized and can be
        % handled by the model, used by the updateState member function.
            varnames = model.primaryVarNames;
            [~, wellVarNames, ~] = ...
                model.FacilityModel.getAllPrimaryVariables(state.wellSol);
            varnames = {varnames{:}, wellVarNames{:}};
        end

        function fds = getAllVarsNames(model)
        % List all the variables that are recognized and can be handled by the model
            fds = {'pressure', 'pressure_matrix', 'wellSol'};
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
        % Updates the state variables that belong to the model, that is the fluid
        % variables. The mechanical variables in the state will be updated
        % by the mechanical model, see DualContMechanicMechModel.
            vars = problem.primaryVariables;
            ind = false(size(vars));
            fluidvars = model.getAllPrimaryVariables(state);

            [lia, loc] = ismember(fluidvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@DualContinuumWaterModel(model, state, problem, ...
                                                                  dx, drivingForces);
        end

    end

end
