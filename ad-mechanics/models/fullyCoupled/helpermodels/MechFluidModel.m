classdef MechFluidModel < ReservoirModel
%
%
% SYNOPSIS:
%   model = MechFluidModel(G, rock, fluid, mech_problem, varargin)
%
% DESCRIPTION:
%   Base class model to set up fully coupled mechanical-fluid simulations.
%   This class is derived for each particular fluid model that is used,
%   see MechBlackOilModel, MechOilWaterModel, MechWaterModel.
%
% PARAMETERS:
%   G            - grid structure
%   rock         - rock structure
%   fluid        - fluid structure
%   mech_problem - Structure that contains the mechanical parameters of the system
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   MechBlackOilModel, MechOilWaterModel, MechWaterModel

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

        % Mechanical model
        mechModel;
        % List of primary variable names for the mechanical part
        MechPrimaryVars;
        % List of all the variable names for the mechanical part
        mechfds;
        % Fluid model
        fluidModel;
        % List of primary variable names for the fluid part
        FluidPrimaryVars;
        % List of all the variable names for the fluid part
        fluidfds;
    end

    methods
        function model = MechFluidModel(G, rock, fluid, mech_problem, varargin)
            model       = model@ReservoirModel(G, varargin{:});
            model.rock  = rock;
            model.fluid = fluid;

            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            % Different fluid models may be used. This base class should be
            % derived for each of those. See e.g. WaterFixedStressFluidModel.m
            % (fixed stress splitting with water phase).
            model.fluidModel = setupFluidModel(model);

            model.fluidfds = model.fluidModel.getAllVarsNames();



            model.mechModel = MechanicMechModel(model.G, rock, mech_problem);
            model.mechfds = model.mechModel.getAllVarsNames();

        end

        function fluidModel = setupFluidModel(model)
            error('Base class function not meant for direct use.');
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)

            error('Base class function not meant for direct use.');
        end

        function [fn, index] = getVariableField(model, name, varargin)
            if ismember(name, model.fluidfds)
                [fn, index] = model.fluidModel.getVariableField(name);
            elseif ismember(name, model.mechfds)
                [fn, index] = model.mechModel.getVariableField(name);
            else
                [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
            end
        end

        function [mechTerm, fluidp] = computeCouplingTerms(model, p0, ...
                                                           xd0, p, xd)
            error('Base class function not meant for direct use.');
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            fluidModel = model.fluidModel;
            mechModel  = model.mechModel;
            [state, fluidReport] = fluidModel.updateState(state, problem, dx, drivingForces);
            [state, mechReport]  = mechModel.updateState(state, problem, dx, []);
            report = [];
        end

        function ctrl = validateDrivingForces(model, ctrl, index)
            fluidmodel = model.fluidModel;
            ctrl = fluidmodel.validateDrivingForces(ctrl, index);
        end

        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel) || isempty(model.FacilityModel.ReservoirModel)
                model.FacilityModel = FacilityModel(model);
                %warning('The MechFluidModel has an empty FacilityModel');
            end
            model.fluidModel.FacilityModel = model.FacilityModel;
            model.fluidModel = model.fluidModel.validateModel(varargin{:});
        end
        
        function [model, state] = updateForChangedControls(model, state, forces)
           [model.fluidModel, state] = model.fluidModel.updateForChangedControls(state, forces);
        end

        function state = validateState(model, state)
           state = model.fluidModel.validateState(state);
           state = model.mechModel.validateState(state);
        end


    end
end
