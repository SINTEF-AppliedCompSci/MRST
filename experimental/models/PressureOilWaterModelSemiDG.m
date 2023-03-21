classdef PressureOilWaterModelSemiDG < PressureOilWaterModel
    
    properties
        disc
    end
    
    methods
        function model = PressureOilWaterModelSemiDG(G, rock, fluid, varargin)
            model= model@PressureOilWaterModel(G, rock, fluid);
            model.disc = [];
            model = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, G.griddim);
            end
        end
        
        % ----------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] ...
                = pressureEquationOilWaterSemiDG(state0, state, model, dt, drivingForces, varargin{:});
        end
        
        % ----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name)
            % Map variables to state field.
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.getVariableField`
            switch(lower(name))
                case {'water', 'swdof'}
                    index = 1;
                    fn = 'sdof';
                case {'oil', 'sodof'}
                    index = 2;
                    fn = 'sdof';
                case {'gas', 'sgdof'}
                    index = 3;
                    fn = 'sdof';
                case{'saturation', 'sdof'}
                    index = ':';
                    fn = 'sdof';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@PressureOilWaterModel(model, name);
            end
        end
        
        % ----------------------------------------------------------------%
        function vars = getSaturationVarNames(model)
            vars = {'sWdof', 'sOdof', 'sGdof'};
            ph = model.getActivePhases();
            vars = vars(ph);
        end
        
        % --------------------------------------------------------------------%
        function state = validateState(model, state)
            % Validate initial state. 
            %
            % SEE ALSO:
            %   :meth:`ad_core.models.PhysicalModel.validateState`
            state = validateState@PhysicalModel(model, state);
            active = model.getActivePhases();
            nPh = nnz(active);
            nc = model.G.cells.num;
            model.checkProperty(state, 'Pressure', [nc, 1], [1, 2]);
%             if nPh > 1
%                 model.checkProperty(state, 'Saturation', [nc, nPh], [1, 2]);
%             end
            if ~isempty(model.FacilityModel)
                state = model.FacilityModel.validateState(state);
            end
        end
    end
    
end

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
