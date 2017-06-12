classdef MechFluidSplitModel < ReservoirModel
% Base class for mechanics-flow splitting. Possibility to implement different
% splitting. For the moment, only fixed-stress splitting is implemented.
%
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

    end

    methods
        function model = MechFluidSplitModel(G, rock, fluid, mech_problem, varargin)

            opt = struct('fluidModelType', 'water', ...
                         'splittingTolerance', 1e-6);
            [opt, rest] = merge_options(opt, varargin{:});
            fluidModelType = opt.fluidModelType;
            model = model@ReservoirModel(G, rest{:});

            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            % Different fluid models may be used. This base class should be
            % derived for each of those. See e.g. WaterFixedStressFluidModel.m
            % (fixed stress splitting with water phase).
            model.fluidModel = model.setupFluidModel(rock, fluid, opt.fluidModelType, ...
                                                           'extraWellSolOutput', ...
                                                           false, rest{:});
            model.fluidfds = model.fluidModel.getAllVarsNames();
            
            
            model.mechModel = MechanicModel(model.G, rock, mech_problem, ...
                                            rest{:});
            model.mechModel.splittingTolerance = opt.splittingTolerance;

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

        function [fn, index] = getVariableField(model, name)
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
                [fn, index] = model.fluidModel.getVariableField(name);
            end
        end

        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel)
                error('The MechFluidSplitModel requires to have an iniatilized FacilityModel')
            end
            model.fluidModel.FacilityModel        = model.FacilityModel;
            
            model = validateModel@ReservoirModel(model, varargin{:});
            return
        end

        
        function state = validateState(model, state)
           state = model.fluidModel.validateState(state);
           state = model.mechModel.validateState(state);
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
