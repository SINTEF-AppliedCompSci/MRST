classdef FluxAddModel < PhysicalModel
    properties
        parentModel
    end
    methods
        function model = FluxAddModel(parent)
            assert(isa(parent, 'PhysicalModel'));
            model = model@PhysicalModel(parent.G);
            model.parentModel = parent;
        end
        function [state, report] = stepFunction(model, state, state0, dt,...
                                                drivingForces, linsolve, nonlinsolve,...
                                                iteration, varargin)
            %state0 = updateConnectionDP(wellmodel, model, state0)
            %state = updateConnectionDP(wellmodel, model, state)
            
            model.parentModel.FacilityModel = FacilityModel(model.parentModel);
            tmpforces = drivingForces;
            tmpforces.W = [];
            wellSol=state.wellSol;
            state.wellSol=[];
            model.parentModel.extraStateOutput=true;
            [problem, state] = model.parentModel.getEquations(state0, state, dt, tmpforces, 'resOnly', true, 'iteration', inf);            
            
            state.wellSol=wellSol;
            for i=1:numel(wellSol)
                %assert(strcmp(drivingForces.W(i).name,wellSol(i).name))
                if( wellSol(i).sign < 0)
                    %wc = vertcat(drivingForces.W(i).cells);
                    wc = wellSol(i).cells;
                    mob = state.mob(wc,:);
                    fm = bsxfun(@rdivide,mob,sum(mob,2));
                    wellSol(i).flux = bsxfun(@times,state.wellSol(i).flux,fm);
                else
                    if(state.wellSol(i).status==1)
                        assert(nnz(state.wellSol(i).compi)==1)
                        wellSol(i).flux=bsxfun(@times,state.wellSol(i).flux,wellSol(i).compi);
                    end
                end
            end
            state.wellSol = wellSol;
            %%{
            %[convergence, values, names] = model.parentModel.checkConvergence(problem);
            convergence = false; values =[]; names = [];
            report = model.makeStepReport(...
                                    'Failure',         false, ...
                                    'Converged',       true, ...
                                    'ResidualsConverged', true(size(convergence)), ...
                                    'Residuals',       values ...
                                    );
            state.convergenceStatus = struct('Residuals', values, 'ResidualsConverged', convergence,...
                                             'Converged', all(convergence), 'Names', {names});
            %}
        end
        
        function [model, state] = updateForChangedControls(model, state, forces)
            % Called whenever controls change. Since this model can be used
            % with wells, we call the facility model's setup routine.
            %model.FacilityModel = model.FacilityModel.setupWells(forces.W);
            %[model.parentModel, state] = model.parentModel.updateForChangedControls(state, forces);
            state.wellSol = initWellSolAD(forces.W, model.parentModel, state);
            [model, state] = updateForChangedControls@PhysicalModel(model, state, forces);
        end
                                    
        function varargout = getActivePhases(model)
            varargout = cell(1, nargout);
            [varargout{:}] = model.parentModel.getActivePhases();
        end
        function forces = getValidDrivingForces(model)
            forces = model.parentModel.getValidDrivingForces();
        end
        function state = validateState(model, state)
            state = model.parentModel.validateState(state);
        end

        

        function model = validateModel(model, varargin)
            %model.parentModel = model.parentModel.validateModel(varargin{:});
            return
        end

        function [fn, index] = getVariableField(model, name)
            [fn, index] = model.parentModel.getVariableField(name);
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
