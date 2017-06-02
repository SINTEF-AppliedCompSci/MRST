classdef MechFluidModel < ReservoirModel
    properties

        % Mechanical model used in the splitting
        mechModel;
        % List of primary variable names for the mechanical part
        MechPrimaryVars;
        % List of all the variable names for the mechanical part
        mechfds;
        % Fluid model used in the splitting
        fluidModel;
        % List of primary variable names for the fluid part
        FluidPrimaryVars;
        % List of all the variable names for the fluid part
        fluidfds;

        alpha_scaling;
        S;
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
            % derived for each of those. See e.g. SinglephaseFixedStressFluidModel.m
            % (fixed stress splitting with single water phase).
            model.fluidModel = setupFluidModel(model);

            model.fluidfds = model.fluidModel.getAllVarsNames();



            model.mechModel = MechanicMechModel(model.G, rock, mech_problem);
            model.mechfds = model.mechModel.getAllVarsNames();

            model.mechModel.alpha_scaling = 1;
            model.mechModel.S = [];
            model.mechModel.ilu_tol = 1e-4;
        end

        function fluidModel = setupFluidModel(model)
            error('Base class function not meant for direct use.');
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)

            error('Base class function not meant for direct use.');
        end

        function [fn, index] = getVariableField(model, name)
            if ismember(name, model.fluidfds)
                [fn, index] = model.fluidModel.getVariableField(name);
            elseif ismember(name, model.mechfds)
                [fn, index] = model.mechModel.getVariableField(name);
            else
                [fn, index] = getVariableField@ReservoirModel(model, name);
            end
        end

        function [mechTerm, fluidp] = computeCouplingTerms(model, p0, ...
                                                           xd0, p, xd)
            error('Base class function not meant for direct use.');
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            fluidModel = model.fluidModel;
            mechModel  = model.mechModel;
            [state, fluidReport] = fluidModel.updateState(state, problem, dx, []);
            [state, mechReport]  = mechModel.updateState(state, problem, dx, []);
            report = []; 
        end

        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel)
                error('The MechFluidModel requires to have an iniatilized FacilityModel')
            end
            model.fluidModel.FacilityModel = model.FacilityModel;
            model = validateModel@ReservoirModel(model, varargin{:});
            return
        end
        
        function state = validateState(model, state)
           state = model.fluidModel.validateState(state);
           state = model.mechModel.validateState(state);
        end
        
        function model = setupOperators(model, G, rock, varargin)


            % Set up divergence / gradient / transmissibility operators for flow
            model = setupOperators@ReservoirModel(model, G, rock, varargin{:});

            operators = setupOperatorsVEM(model.G, model.mech.el_bc, ...
                                                   model.mech.load, ...
                                                   model.alpha_scaling, ...
                                                   model.S);
            model.operators.mech = operators.mech;
            model.operators.extra = operators.extra;

        end

    end
end
