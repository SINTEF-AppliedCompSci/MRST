classdef BiotModel < PhysicalModel

    properties
        
        rock
        mech
        fluid
        
        % Property container
        MechBiotPropertyFunctions
        FluidBiotPropertyFunctions 
    end

    methods
        
        function model = BiotModel(G, rock, fluid, mech, varargin)
            
            model = model@PhysicalModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
            
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end
            % Physical properties of rock and fluid
            model.mech  = mech;
            model.fluid = fluid;
            
            % Add mechanical operators  
            model.operators = setupBiotOperators(model);

            % Add mechanical properties
            model.MechBiotPropertyFunctions = MechBiotPropertyFunctions(model);
            model.FluidBiotPropertyFunctions = FluidBiotPropertyFunctions(model);

        end
        

        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@PhysicalModel(model);
            containers = [containers, {model.MechBiotPropertyFunctions, model.FluidBiotPropertyFunctions}];
        end

        function model = setupStateFunctionGroupings(model, varargin) 
            
            model = setupStateFunctionGroupings@PhysicalModel(model, varargin{:});
            
            % Use mechanical coupling in flow equations
            fp = model.FlowPropertyFunctions;
            fp = fp.setStateFunction('BasePoreVolume', BlackOilPoreVolume(model));
            fp = fp.setStateFunction('PoreVolume'    , BiotPoreVolume(model));
            fp = fp.setStateFunction('Dilatation'    , BiotCoupledDilatation(model));
            model.FlowPropertyFunctions = fp;
            
            % Use fluid coupling in mechanical equations
            mp = model.MechBiotPropertyFunctions;
            mp = mp.setStateFunction('BiotGradP', BiotCoupledGradP(model));
            model.MechBiotPropertyFunctions = mp;
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            [eqs, names, types, state] = biotEquations(model, state0, state, dt, drivingForces, varargin)
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)

            [vars, names, origin] = getPrimaryVariables@PhysicalModel(model, state);
            [u, lm, p, lf] = model.getProps(state, 'u', 'lambdamech', 'pressure', 'lambdafluid');
            vars = [vars, {u, lm, p, lf}];
            names = [names, {'displacement', 'lambdamech', 'pressure', 'lambdafluid'}];

            origin = [origin, {class(model)}];
        end
        
        function state = validateState(model, state)

            state = validateState@PhysicalModel(model, state);
            if ~isfield(state, 'wellSol')
                state.wellSol = []; % (needed by simulateScheduleAD)
            end
            
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            % Support for wells (needed by simulateScheduleAD)
            forces.W   = [];
        end
        
        function [fn, index] = getVariableField(model, name, varargin)
            
            switch(lower(name))
              case {'u', 'displacement'}
                % Displacement
                fn = 'u';
                index = ':';
              case {'lambdamech'}
                % Lagrangian variables for corresponding to the boundary conditions
                fn = lower(name);
                index = ':';
              case {'pressure', 'p'}
                fn = 'pressure';
                index = ':';
              case {'lambdafluid'}
                fn = lower(name);
                index = ':';
              otherwise
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end

        end

    end

end
