classdef BiotMechModel < PhysicalModel

    properties
        

        mech;
        % Structure that contains the mechanical parameter of the system
        % This structure should contain the fields prop
        % prop :  - lambda, first Lame parameter
        %         - mu    , second Lame parameter
        %
        % and the loading structure loadstruct
        % loadstruct:  - bc      , Dirichlet boundary conditions 
        %              - extforce, Neumann boundary conditions (external forces)
        %              - force   , Volumetric boundary conditions
        
        
        rock;
        % Structure that contains the rock properties. The structure must
        % contain the Biot coefficient
        %
        % alpha - Biot coefficient (one entry per cell)
        %
        
        % Property container
        MechBiotPropertyFunctions;
        
    end

    methods
        function model = BiotMechModel(G, rock, mech, varargin)
            
            model = model@PhysicalModel(G, varargin{:});

            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            % Physical properties of rock and fluid
            model.mech  = mech;
            model.rock  = rock;

            model.MechBiotPropertyFunctions = MechBiotPropertyFunctions(model);

            model.operators = setupMpsaOperators(model);
            
        end

        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            eqs{1} = model.getProp(state, 'MomentumEquations');
            names = {'momentum'};
            types = {'mixed'};
            
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)
            [u, lambdamech] = model.getProps(state, 'displacement', 'lambdamech');
            vars = {u, lambdamech};
            names = {'displacement', 'lambdamech'};
            origin = {class(model)};
        end
        
        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@PhysicalModel(model);
            containers = [containers, {model.MechBiotPropertyFunctions}];
        end
        
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@PhysicalModel(model, varargin{:});
            model.MechBiotPropertyFunctions = MechBiotPropertyFunctions(model);
        end
    
        function [fn, index] = getVariableField(model, name, varargin)
        % Get the index/name mapping for the model (such as where pressure or water saturation is located in state)
            switch(lower(name))
              case {'u', 'displacement'}
                % Displacement
                fn = 'u';
                index = ':';
              case {'lambdamech'}
                % Lagrangian variables for corresponding to the boundary conditions
                fn = 'lambdamech';
                index = ':';
              case {'biotgradp'}
                % external force corresponding to a pressure gradient.
                fn = 'biotgradp';
                index = ':';
              otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end
        end

    end
end
