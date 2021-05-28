classdef TracerModel < ReservoirModel
    properties
        tracerNames
    end

    methods
        function model = TracerModel(G, rock, fluid, varargin)
            model = model@ReservoirModel(G, rock, fluid);
            model = merge_options(model, varargin{:});
            model.stepFunctionIsLinear = true;
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        % see function equationsTracer for specification of drivingForces.bc
            [problem, state] = equationsTracer(state0, state, model, dt, ...
                                               drivingForces, varargin{:});

        end

        function n = getNumberOfTracers(model)
            n = numel(model.tracerNames);
        end
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@PhysicalModel(model, state0, state, dt, drivingForces);
            if ~isempty(model.FacilityModel)
                %state.wellSol = model.FacilityModel.updateWellSolAfterStep(state.wellSol, state0.wellSol);
            end
        end
        function [model, state] = updateForChangedControls(model, state, forces)
            % Called whenever controls change. Since this model can be used
            % with wells, we call the facility model's setup routine.
            %model.FacilityModel = model.FacilityModel.setupWells(forces.W);
            %[model.parentModel, state] = model.parentModel.updateForChangedControls(state, forces);
            state.wellSol = initWellSolAD(forces.W, model, state);
            [model, state] = updateForChangedControls@PhysicalModel(model, state, forces);
        end
        
        function state = validateState(model, state)
            % Check parent class
            %state = validateState@ReservoirModel(model, state);
            nc = model.G.cells.num;
            if isfield(state, 'tracer')
                model.checkProperty(state, 'tracer', nc, 1);
            else
                nt = model.getNumberOfTracers();
                state.tracer = zeros(nc, nt);
            end
        end
        
        function model = validateModel(model, varargin)
            if isempty(model.FacilityModel)
                model.FacilityModel = FacilityModel(model);
            end
            if nargin > 1
                model.FacilityModel = model.FacilityModel.setReservoirModel(model);
                W = varargin{1}.W;
                model.FacilityModel.WellModels = {};
                model.FacilityModel = model.FacilityModel.setupWells(W);
            end
            model = validateModel@PhysicalModel(model, varargin{:});
            return
        end
   
        
        
        function [fn, index] = getVariableField(model, name, varargin)
            tsub = strcmpi(model.tracerNames, name);
            if any(tsub)
                fn = 'tracer';
                index = find(tsub);
            else
                switch(lower(name))
                    case 'tracer'
                        fn = 'tracer';
                        index = ':';
                    otherwise
                        % Basic phases are known to the base class
                        [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
                end
            end
        end
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
