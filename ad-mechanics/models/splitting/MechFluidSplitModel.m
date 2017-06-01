classdef MechFluidSplitModel < ReservoirModel
% Base class for mechanics-flow splitting.
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

        %
        alpha_scaling;
        S;
        ilu_tol;
    end

    methods
        function model = MechFluidSplitModel(G, rock, fluid, mech_problem, varargin)

            opt = struct('fluidModelType', 'single phase');
            [opt, rest] = merge_options(opt, varargin{:});
            fluidModelType = opt.fluidModelType;

            model = model@ReservoirModel(G, rest{:});

            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            % Different fluid models may be used. This base class should be
            % derived for each of those. See e.g. SinglephaseFixedStressFluidModel.m
            % (fixed stress splitting with single water phase).
            model.fluidModel = setupFluidModel(model, rock, fluid, ...
                                                      opt.fluidModelType, ...
                                                      'extraWellSolOutput', ...
                                                      false, rest{:});

            model.fluidfds = model.fluidModel.getListFields();

            model.alpha_scaling = 1;
            model.S = [];
            model.ilu_tol = 1e-4;
            model = merge_options(model, rest{:});

            model.mechModel = MechanicModel(model.G, rock, mech_problem);
            model.mechfds = model.mechModel.getListFields();

            model.mech_solver = NonLinearSolver();
            model.fluid_solver = NonLinearSolver();


            %     switch model.fluidModelType
            %       case 'single phase'
            %         model.water = true;
            %         model.oil = false;
            %         model.gas = false;
            %         model.saturationVarNames = {};
            %         model.fluidfds = {'wellSol', 'pressure'};
            %         model.mechfds = {'xd', 'uu', 'u', 'stress', 'strain', 'vdiv'};
            %       case 'oil water'
            %         model.water = true;
            %         model.oil = true;
            %         model.gas = false;
            %         model.saturationVarNames = {'sw', 'so'};
            %         model.fluidfds = {'wellSol', 'pressure', 's'};
            %         model.mechfds = {'xd', 'uu', 'u', 'stress', 'strain', 'vdiv'};
            %       case 'blackoil'
            %         model.water = true;
            %         model.oil = true;
            %         model.gas = true;
            %         model.disgas = true;
            %         model.vapoil = false;
            %         model.saturationVarNames = {'sw', 'so', 'sg'};
            %         model.fluidfds = {'wellSol', 'pressure', 's', 'rs', 'rv'};
            %         model.mechfds = {'xd', 'uu', 'u', 'stress', 'strain', 'vdiv'};
            %       otherwise
            %         error('fluidModelType not recognized.')
            %     end

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
            model.fluidModel.FacilityModel = model.FacilityModel;
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
