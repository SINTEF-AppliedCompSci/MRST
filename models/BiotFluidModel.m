classdef BiotFluidModel < PhysicalModel

    properties
         
        rock
        fluid
        
        FluidBiotPropertyFunctions % Grouping for flow properties
    end
    
    methods
        
        function model = BiotFluidModel(G, rock, fluid, varargin)
            
            model = model@PhysicalModel(G, varargin{:});
            model = merge_options(model, varargin{:});
            
            model.rock     = rock;
            model.fluid    = fluid;
            
            model.operators = setupMpfaOperators(model);
            model.FluidBiotPropertyFunctions = FluidBiotPropertyFunctions(model);            
        end
        
        function [vars, names, origin] = getPrimaryVariables(model, state)
            [p, lambdafluid] = model.getProps(state, 'pressure', 'lambdafluid');
            vars = {p, lambdafluid};
            names = {'pressure', 'lambdafluid'};
            origin = {class(model)};
        end
        
        
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@PhysicalModel(model);
            containers = [containers, {model.FluidBiotPropertyFunctions}];
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            
            model = setupStateFunctionGroupings@PhysicalModel(model, varargin{:});
            
            fp = model.FluidBiotPropertyFunctions;
            % fp = fp.setStateFunction('BasePoreVolume', BlackOilPoreVolume(model));
            fp = fp.setStateFunction('PoreVolume', PoreVolume(model));
            % fp = fp.setStateFunction('Dilatation'    , BiotDilatation(model));
            
            model.FluidBiotPropertyFunctions = fp;
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            
            [eqs, names, types, states] = fluidEquations(model, state0, state, dt, drivingForces, varargin{:});
            
        end

        function state = validateState(model, state)
            state = validateState@PhysicalModel(model, state);
            if ~isfield(state, 'wellSol')
                state.wellSol = [];
            end
        end
                
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            % Support for wells
            forces.W   = [];
        end


        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
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

