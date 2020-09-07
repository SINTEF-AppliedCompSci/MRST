classdef DualContMechFluidModel < DualContinuumReservoirModel
%
% SYNOPSIS:
%   model = DualContMechFluidModel(G, rock, fluid, rock_matrix, fluid_matrix, 
%                                  mech_problem, varargin)
%
% DESCRIPTION: 
%   Base class model to set up fully coupled dual-permeability mechanical-fluid
%   simulations. This class is derived for each particular fluid model that is
%   used, see DualContMechWaterModel.
%
% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock cell for fracture / matrix
%   fluid        - Fluid cell for fracture / matrix
%   mech_problem - Structure that contains the mechanical parameters of the system
%
% RETURNS:
%   class instance
%
% EXAMPLE: 
%
% SEE ALSO: DualContMechWaterModel
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
        % Grouping for flow properties
        DC_PoroelasticPropertyFunctions; 
    end

    methods
        function model = DualContMechFluidModel(G, rock, fluid, mech_problem, varargin)
            model = model@DualContinuumReservoirModel(G, rock, fluid, varargin{:});
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            % Different fluid models may be used. 
            model.fluidModel = setupFluidModel(model); 
            model.fluidModel.extraStateOutput = 'True';
            model.fluidfds = model.fluidModel.getAllVarsNames(); % get variable names for fluid proble

            % Setup mech model
            rock_fracture = rock{1};
            rock_matrix = rock{2};
            model.mechModel = DualContMechMechanicModel(model.G, rock_fracture, rock_matrix, mech_problem); % Mechanical problem for FC sim
            model.mechfds = model.mechModel.getAllVarsNames(); % get variable names for mechanical problem
            model.DC_PoroelasticPropertyFunctions = [];
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
                [fn, index] = getVariableField@DualContinuumReservoirModel(model, name, varargin{:});
            end
        end

        function mechTerm = computeStrainTerms(model, xd0, xd)
            error('Base class function not meant for direct use.');
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            fluidModel = model.fluidModel; 
            mechModel  = model.mechModel;
            [state, fluidReport] = fluidModel.updateState(state, problem, dx, []);
            [state, mechReport]  = mechModel.updateState(state, problem, dx, []);
            report = [];
            state = addDerivedQuantitiesDC(model, state);
        end

        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel)
                warning('The DualContMechFluidModel has an empty FacilityModel');
            end
            % let the fluid model deal with the FacilityModel
            model.fluidModel.FacilityModel = model.FacilityModel;
            model.fluidModel = model.fluidModel.validateModel(); % setup wells 
            if isempty(model.DC_PoroelasticPropertyFunctions)
                model.DC_PoroelasticPropertyFunctions = DC_PoroelasticPropertyFunctions(model); %#ok
            end
            return
        end

        function state = validateState(model, state)
           state = model.fluidModel.validateState(state);
           state = model.mechModel.validateState(state);
        end

        function containers = getStateFunctionGroupings(model)
            % called in getEquations when calling
            % model.initStateFunctionContainers() method
            containers = getStateFunctionGroupings@PhysicalModel(model);
            extra = {model.DC_PoroelasticPropertyFunctions};
            if ~isempty(model.FacilityModel)
                fm_props = model.FacilityModel.getStateFunctionGroupings();
                extra = [extra, fm_props];
            end
            extra = extra(~cellfun(@isempty, extra));
            containers = [containers, extra];
        end
    end
end
