classdef MandelModel < BiotModel


    properties
        topfaces
    end
    

    methods
        
        function model = MandelModel(G, rock, fluid, mech, topfaces, varargin)
            
            model = model@BiotModel(G, rock, fluid, mech, varargin{:});
            model = merge_options(model, varargin{:});
            model.topfaces = topfaces;
            model.operators = setMandelOperators(model);
            
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
        
        function operators = setMandelOperators(model)
            
            op = model.operators;
            topfaces = model.topfaces;
            tbls = model.operators.tbls;
            bc = model.mech.loadstruct.bc;
            
            % We recover the degrees of freedom that belongs to the top surface
            bcnodefacetbl = bc.bcnodefacetbl;
            bcnodefacetbl = bcnodefacetbl.addLocInd('bcinds');

            topfacetbl.faces = topfaces;
            topfacetbl = IndexArray(topfacetbl);

            map = TensorMap();
            map.fromTbl = topfacetbl;
            map.toTbl = bcnodefacetbl;
            map.mergefds = {'faces'};
            map = map.setup();

            toplambdamech = map.eval(ones(topfacetbl.num, 1));
            toplambdamech = find(toplambdamech);
            
            
            
            
        end
        

    end

end
