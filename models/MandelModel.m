classdef MandelModel < BiotModel



    methods
        
        function model = MandelModel(G, rock, fluid, mech, varargin)
            
            model = model@BiotModel(G, rock, fluid, mech, varargin{:});
            model = merge_options(model, varargin{:});
            
        end
        

        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces, varargin)
            [eqs, names, types, state] = mandelEquations(model, state0, state, dt, drivingForces, varargin);
        end

        function [vars, names, origin] = getPrimaryVariables(model, state)

            [vars, names, origin] = getPrimaryVariables@PhysicalModel(model, state);
            [u, p, lm, lf, vd] = model.getProps(state, 'u', 'pressure', 'lambdamech', 'lambdafluid', 'vertdisp');
            vars = [vars, {u, p, lm, lf, vd}];
            names = [names, {'displacement', 'pressure', 'lambdamech', 'lambdafluid', 'vertdisp'}];

            origin = [origin, {class(model)}];
        end
        

        function [model, state] = prepareReportstep(model, state, state0, dT, drivingForces)
            [model, state] = prepareReportstep@PhysicalModel(model, state, state0, dT, drivingForces);
            extforce = drivingForces.extforce;
            state = model.setProp(state, 'avgtopforce', avgtopforce);
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            % Support for wells (needed by simulateScheduleAD)
            forces.W   = [];
            % mechanical forces
            forces.avgtopforce = []; 
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
              case {'avgtopforce', 'F'}
                % Average of the vertical component of the force exterted at the top.
                fn = 'F';
                index = ':';
              case {'vertdisp', 'vd'}
                % vertical displacement at the top is constant and equal to scalar stored in 'vd'
                fn = 'vd';
                index = ':';
              otherwise
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
            end

        end

    end

end
