classdef BiotFluidModel < PhysicalModel

    properties
        rock
        FluidBiotPropertyFunctions % Grouping for flow properties
    end
    
    methods
        function model = BiotFluidModel(G, rock, fluid, varargin)
            model = model@PhysicalModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
            
            operators = setupMpfaOperators(model);
            
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@PhysicalModel(model, varargin{:});
            
            fp = model.FlowPropertyFunctions;
            fp = fp.setStateFunction('BasePoreVolume', BlackOilPoreVolume(model));
            fp = fp.setStateFunction('PoreVolume'    , BiotPoreVolume(model));
            fp = fp.setStateFunction('Dilatation'    , BiotDilatation(model));
            
            model.FlowPropertyFunctions = fp;
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            
            eqs = fluidequations(model, state0, state, dt, drivingForces, varargin{:});
            names = {'fluid'};
            types = {'mixed'};
            
        end
 
        function [fn, index] = getVariableField(model, name, varargin)
            switch(lower(name))
              case {'pressure', 'p'}
                fn = lower(name);
                index = ':';
              case {'lambda'}
                fn = lower(name);
                index = ':';
              otherwise
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end
            
        end

    end
    
end

