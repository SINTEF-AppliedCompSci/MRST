classdef MechFluidSplitModel < ReservoirModel
%
%
% SYNOPSIS:
%   model = MechFluidSplitModel(G, rock, fluid, mech_problem, varargin)
%
% DESCRIPTION:
%   Base class for mechanics-flow splitting models. It offers the
%   possibility to implement different splitting. For the moment, only
%   fixed-stress splitting is implemented. The model contains a mechanical
%   and fluid model, which can independently handle the mechanical and
%   fluid systems of equations.
%
% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock structure
%   fluid        - Fluid structure
%   mech_problem - Structure that contains the mechanical parameters of the system
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   MechFluidFixedStressSplitModel

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
        % Mechanical model used in the splitting
        mechModel;
        % List of variable names for the mechanical part (see sync functions below)
        mechfds;
        % Fluid model used in the splitting
        fluidModel;
        % List of variable names for the fluid part (see sync functions below)
        fluidfds;
        % Solver to be used for the mechanical part
        mech_solver;
        % Solver to be used for the fluid part
        fluid_solver;

        % Tolerance used in the splitting scheme
        splittingTolerance
        % Splitting verbose
        splittingVerbose
    end

    methods
        function model = MechFluidSplitModel(G, rock, fluid, mech_problem, varargin)

            opt = struct('fluidModelType', 'water', ...
                         'splittingTolerance', 1e-6, ...
                         'splittingVerbose', false);
            [opt, rest] = merge_options(opt, varargin{:});
            fluidModelType = opt.fluidModelType;

            model = model@ReservoirModel(G, rest{:});
            if opt.splittingVerbose
                model.verbose = true;
            end

            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            model.water = ...
               ~isempty(regexp(fluidModelType, '\<(water|blackoil)\>', 'once'));

            model.oil = ...
               ~isempty(regexp(fluidModelType, '\<(oil|blackoil)\>', 'once'));

            model.gas = ...
               ~isempty(regexp(fluidModelType, '\<blackoil\>', 'once'));

            % Different fluid models may be used. This base class should be
            % derived for each of those. See e.g. WaterFixedStressFluidModel.m
            % (fixed stress splitting with water phase).
            model.fluidModel = model.setupFluidModel(rock, fluid, ...
                                                     fluidModelType, ...
                                                     'extraWellSolOutput', ...
                                                     false, rest{:});
            model.fluidfds = model.fluidModel.getAllVarsNames();

            model.mechModel = MechanicModel(model.G, rock, mech_problem, ...
                                            rest{:});
            model.splittingTolerance = opt.splittingTolerance;

            model.mechfds = model.mechModel.getAllVarsNames();

            model.mech_solver = NonLinearSolver();
            model.fluid_solver = NonLinearSolver();
        end

        function fluidModel = setupFluidModel(model, rock, fluid, fluidModelType, ...
                                                     varargin)
            error('Base class function not meant for direct use.');
        end

        function [state, report] = stepFunction(model, state, state0, dt, ...
                                                drivingForces, linsolve, ...
                                                nonlinsolve, iteration, ...
                                                varargin)
            error('Base class function not meant for direct use.');
        end

        function divTerm = computeMechTerm(model, state)
            error('Base class function not meant for direct use.');
        end

        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
              case {'xd'}
                fn = 'xd';
                index = 1;
              case {'uu'}
                fn = 'uu';
                index = ':';
              case {'u'}
                fn = 'u';
                index = ':';
              case {'stress'}
                fn = 'stress';
                index = ':';
              case {'strain'}
                fn = 'strain';
                index = ':';
              case {'vdiv'}
                fn = 'vdiv';
                index = ':';
              otherwise
                % This will throw an error for us
                [fn, index] = model.fluidModel.getVariableField(name, varargin{:});
            end
        end

        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel) || isempty(model.FacilityModel.ReservoirModel)
               model.FacilityModel = FacilityModel(model);
               %error('The MechFluidSplitModel requires to have an iniatilized FacilityModel')
            end
            model.fluidModel.FacilityModel        = model.FacilityModel;
            if isempty(model.operators)
                model.operators = model.fluidModel.operators;
            end
            model = validateModel@ReservoirModel(model, varargin{:});
            return
        end


        function [convergence, values, names] = checkConvergence(model, problem, ...
                                                              varargin)

            mechModel = model.mechModel;
            mechModel.nonlinearTolerance = model.splittingTolerance;

            [convergence, values, names] = mechModel.checkConvergence(problem);
        end


        function state = validateState(model, state)
           state = model.fluidModel.validateState(state);
           state = model.mechModel.validateState(state);
        end


        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, ...
                                                              state0, state, dt, drivingForces);
            state = addDerivedQuantities(model.mechModel, state);
        end

        function wstate = syncWStateFromState(model, state)
            wstate = updateFields(model.fluidModel, [], model, state, model.fluidfds);
        end

        function mstate = syncMStateFromState(model, state)
            mstate = updateFields(model.mechModel, [], model, state, model.mechfds);
        end

        function state = syncStateFromWState(model, state, wstate)
            state = updateFields(model, state, model.fluidModel, wstate, model.fluidfds);
        end

        function state =syncStateFromMState(model, state, mstate)
            state = updateFields(model, state, model.mechModel, mstate, model.mechfds);
        end

    end
end

function outState = updateFields(outModel, outState, inModel, inState, flds)
    for i = 1 : numel(flds)
        val      = inModel.getProp(inState, flds{i});
        outState = outModel.setProp(outState, flds{i}, val);
    end
end
