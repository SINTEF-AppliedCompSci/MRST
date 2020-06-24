classdef BiotBlackOilModel < GenericBlackOilModel

    properties
        
        mech

        eta
        bcetazero
        
        BiotPropertyFunctions
        
    end

    methods
        
        function model = BiotBlackOilModel(G, rock, fluid, mech, varargin)
            
            model = model@GenericBlackOilModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
            
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end
            % Physical properties of rock and fluid
            model.mech  = mech;
        
            model.eta = 1/3;
            model.bcetazero = false;
            
            % Add mechanical operators  
            model.operators = setupBiotOperators(model);

            % Add mechanical properties
            model.MechBiotPropertyFunctions = MechBiotPropertyFunctions(model);
            model.FluidBiotPropertyFunctions = FluidBiotPropertyFunctions(model);

        end
        

        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@PhysicalModel(model);
            containers = [containers, {model.BiotPropertyFunctions}];
        end

        function model = setupStateFunctionGroupings(model, varargin) 
            
            model = setupStateFunctionGroupings@GenericBlackOilModel(model, varargin{:});
            
            biotprops = model.BiotPropertyFunctions; 
            pvtprops  = model.PVTPropertyFunctions; 
            
            pv = pvtprops.getStateFunction('PoreVolume');
            biotprops = biotprops.setStateFunction('BasePoreVolume', pv);
            pvtprops  = pvtprops.setStateFunction('PoreVolume', BiotPoreVolume(model));
            biotprops = biotprops.setStateFunction('Dilatation', BiotCoupledDilatation(model));
            
            model.BiotPropertyFunctions = biotprops;
            model.PVTPropertyFunctions  = pvtprops;
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            [eqs, names, types, state] = biotBlackOilEquations(model, state0, state, dt, drivingForces, varargin);
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)

            [vars, names, origin] = getPrimaryVariables@GenericBlackOilModel(model, state);
            [u, bp, lm, lf] = model.getProps(state, 'u', 'biotpressure', 'lambdamech', 'lambdafluid');
            vars = [vars, {u, p, lm, lf}];
            names = [names, {'displacement', 'biotpressure', 'lambdamech', 'lambdafluid'}];

            origin = [origin, {class(model)}];
        end
        
        function state = validateState(model, state)

            state = validateState@PhysicalModel(model, state);
            if ~isfield(state, 'wellSol')
                state.wellSol = []; % (needed by simulateScheduleAD)
            end
            
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
              case {'biotpressure', 'bp'}
                fn = 'biotpressure';
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
