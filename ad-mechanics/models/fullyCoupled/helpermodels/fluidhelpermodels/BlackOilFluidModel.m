classdef BlackOilFluidModel < ThreePhaseBlackOilModel
%
%
% SYNOPSIS:
%   model = BlackOilFluidModel(G, rock, fluid, varargin)
%
% DESCRIPTION:
%   This model is for the fluid part of a fully coupled poroelastic
%   simulation. It is derived from the fluid model
%   ThreePhaseBlackOilModel and some few functionalities that are
%   needed for the coupled solver are added.
%
% PARAMETERS:
%   G        - Grid structure
%   rock     - Rock structure
%   fluid    - Fluid structure
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   MechFluidModel

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

    properties
        % Primary variables that are used in the fluid system
        primaryVarNames;
    end

    methods
        function model = BlackOilFluidModel(G, rock, fluid, varargin)
        % Constructor
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);

            model.disgas          = true;
            model.vapoil          = false;
            model.primaryVarNames = {'pressure', 'sW', 'x'}; % well variables
                                                             % not included
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            error(['The BlackOilFluidModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from MechFluidModel.'])
        end

        function varnames = getAllPrimaryVariables(model, state)
        % list all the primarty variables that are recognized and can be
        % handled by the model, used by the updateState member function.
            varnames = model.primaryVarNames;
            [~, wellVarNames, ~] = ...
                model.FacilityModel.getAllPrimaryVariables(state.wellSol);
            varnames = {varnames{:}, wellVarNames{:}};
        end

        function fds = getAllVarsNames(model)
        % list all the variables that are recognized and can be handled by the model
            fds = {'wellSol', 'pressure', 's', 'rs', 'rv', 'sW', 'sG', 'sO', ...
                   'water', 'oil', 'gas'};
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            % updates the state variable that belongs to the model, that is the fluid
            % variables. The mechanical variables in the state will be updated
            % by the mechanical model, see MechanicMechModel.
            vars = problem.primaryVariables;
            ind = false(size(vars));
            fluidvars = model.getAllPrimaryVariables(state);

            [lia, loc] = ismember(fluidvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@ThreePhaseBlackOilModel(model, state, ...
                                                              problem, dx, ...
                                                              drivingForces);

        end


    end

end
