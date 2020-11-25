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
            avgtopforce = drivingForces.avgtopforce;
            state = model.setProp(state, 'avgtopforce', avgtopforce);
        end
        
        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@PhysicalModel(model);
            % Support for wells (needed by simulateScheduleAD)
            forces.W   = [];
            % mechanical forces
            forces.avgtopforce = []; 
        end
        
        function model = setupStateFunctionGroupings(model, varargin) 
            
            model = setupStateFunctionGroupings@BiotModel(model, varargin{:});
           
            biotprops = model.BiotPropertyFunctions;
            biotprops = biotprops.setStateFunction('Dilatation', BiotBlackOilDilatation(model));
            model.BiotPropertyFunctions = biotprops;
            
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
            
            operators = model.operators;
            topfaces  = model.topfaces;
            mech      = model.mech;
            G         = model.G;
            
            tbls = operators.tbls;
            nodefacetbl = tbls.nodefacetbl;
            facetbl = tbls.facetbl;
            
            bc = model.mech.loadstruct.bc;
            
            % We recover the degrees of freedom that belong to the top surface
            bcnodefacetbl = bc.bcnodefacetbl;
            bcnodefacetbl = bcnodefacetbl.addLocInd('bcinds');

            topfacetbl.faces = topfaces;
            topfacetbl = IndexArray(topfacetbl);

            bctopnodefacetbl = crossIndexArray(topfacetbl, bcnodefacetbl, {'faces'});
            
            map = TensorMap();
            map.fromTbl = nodefacetbl;
            map.toTbl = facetbl;
            map.mergefds = {'faces'};
            map = map.setup();
            
            nnodeperface = map.eval(ones(nodefacetbl.num, 1));
            nfareas = 1/nnodeperface.*(G.faces.areas);
            
            map = TensorMap();
            map.fromTbl = facetbl;
            map.toTbl = bctopnodefacetbl;
            map.mergefds = {'faces'};
            map = map.setup();
            
            R = map.eval(nfareas);
            
            map = TensorMap();
            map.fromTbl = bctopnodefacetbl;
            map.toTbl = bcnodefacetbl;
            map.mergefds = {'nodes', 'faces', 'bcinds'};
            map = map.setup();
            
            R = map.eval(R);
            
            operators.R = R;
            
            % remove dependency in extforce for the divergence operator
            
            divuopext = operators.divuop;
            extforce = mech.loadstruct.extforce;
            divuop = @(u, p, lm) divuopext(u, p, lm, extforce);
            operators.divuop = divuop;
            
            % same for momentop
            momentopext = operators.momentop;
            extforce = mech.loadstruct.extforce;
            momentop = @(u, p, lm) momentopext(u, p, lm, extforce);
            operators.momentop = momentop;
            
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

