classdef TransportOilWaterModelDG < TransportOilWaterModel
    % Two phase oil/water system without dissolution
    properties
        degree
    end

    methods
        function model = TransportOilWaterModelDG(G, rock, fluid, varargin)
            model = model@TransportOilWaterModel(G, rock, fluid);
            model.degree = 1;
            
            model = merge_options(model, varargin{:});
            
            cellInt = makeCellIntegrator(model.G, (1:G.cells.num)', model.degree, 'tri');
            model.operators.cellInt = cellInt;
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = transportEquationOilWaterDG(state0, state, model,...
                               dt, ...
                               drivingForces,...
                               'solveForOil',   model.conserveOil, ...
                               'solveForWater', model.conserveWater, ...
                               varargin{:});
            
        end
        
        function [fn, index] = getVariableField(model, name)
        % Map variables to state field.
        %
        % SEE ALSO:
        %   :meth:`ad_core.models.PhysicalModel.getVariableField`
        switch(lower(name))
            case {'water'}
                index = ':';
                fn = 'sdof';
            otherwise
                % This will throw an error for us
                [fn, index] = getVariableField@TransportOilWaterModel(model, name);
        end
    end
        
        
        
%         function p = getProp(model, state, name)
%             switch(lower(name))
%                 case {'sw', 'water'}
%                     
%                     [k, nDof] = dgBasis(model.degree, model.G.griddim);
%                     p = @(x) x*0;
%                     for dofNo = 1:nDof
%                         psi = Polynomial(k(dofNo,:));
%                         p  = @(x) p(x) + state.sWdof(:, dofNo)*psi(x);
%                     end  
%                     
%                 otherwise
%                     p = getProp@TransportOilWaterModel(model, state, name);
%             end    
%         end
                    
                    
        
        
        
        
    end
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}