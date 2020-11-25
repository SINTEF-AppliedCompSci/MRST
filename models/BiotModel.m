classdef BiotModel < PhysicalModel

    properties
        rock
        mech
        fluid

        eta
        bcetazero

        % Property container
        BiotPropertyFunctions
        MechPropertyFunctions
    end

    methods
        
        function model = BiotModel(G, rock, fluid, mech, varargin)
            
            model = model@PhysicalModel(G, varargin{:});
            model = merge_options(model, varargin{:});

            % Physical properties of rock and fluid
            model.rock  = rock;
            model.mech  = mech;
            model.fluid = fluid;
        
            model.eta = 1/3;
            model.bcetazero = false;
            
            % setup operators
            model.operators = setupBiotOperators(model);

            % setup properties
            model.BiotPropertyFunctions = BiotPropertyFunctions(model);
            model.MechPropertyFunctions = MechPropertyFunctions(model);
            
        end
        

        function containers = getStateFunctionGroupings(model)
            containers = getStateFunctionGroupings@PhysicalModel(model);
            containers = [containers, {model.BiotPropertyFunctions, model.MechPropertyFunctions}];
        end

        function model = setupStateFunctionGroupings(model, varargin) 
            
            model = setupStateFunctionGroupings@PhysicalModel(model, varargin{:});
            
            bp = model.BiotPropertyFunctions;
            mp = model.MechPropertyFunctions;
            
            bp = bp.setStateFunction('Dilatation', BiotDilatation(model));
            mp = mp.setStateFunction('FaceNodeDisplacement', FaceNodeDisplacementExtforce(model));
            
            model.BiotPropertyFunctions = bp;
            model.MechPropertyFunctions = mp;
            
        end
        
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            [eqs, names, types, state] = biotEquations(model, state0, state, dt, drivingForces, varargin);
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)

            [vars, names, origin] = getPrimaryVariables@PhysicalModel(model, state);
            [u, p, lm, lf] = model.getProps(state, 'u', 'pressure', 'lambdamech', 'lambdafluid');
            vars = [vars, {u, p, lm, lf}];
            names = [names, {'displacement', 'pressure', 'lambdamech', 'lambdafluid'}];

            origin = [origin, {class(model)}];
        end
        
        function state = validateState(model, state)

            state = validateState@PhysicalModel(model, state);
            if ~isfield(state, 'wellSol')
                state.wellSol = []; % (needed by simulateScheduleAD)
            end
            
        end

        function [model, state] = prepareReportstep(model, state, state0, dT, drivingForces)
            [model, state] = prepareReportstep@PhysicalModel(model, state, state0, dT, drivingForces);
            extforce = drivingForces.extforce;
            state = model.setProp(state, 'extforce', extforce);
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            % Support for wells (needed by simulateScheduleAD)
            forces.W   = [];
            % mechanical forces
            forces.extforce = []; 
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
              case {'extforce'}
                % external force
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

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}

