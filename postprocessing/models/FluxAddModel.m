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
            for i=1:numel(drivingForces.W)
                assert(strcmp(drivingForces.W(i).name,wellSol(i).name))
                if( drivingForces.W(i).sign < 0)
                    wc = vertcat(drivingForces.W(i).cells);
                    mob = state.mob(wc,:);
                    fm = bsxfun(@rdivide,mob,sum(mob,2));
                    wellSol(i).flux = bsxfun(@times,state.wellSol(i).flux,fm);
                else
                    assert(nnz(state.wellSol(i).compi)==1)
                    wellSol(i).flux=bsxfun(@times,state.wellSol(i).flux,drivingForces.W(i).compi);
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

        function [model, state] = updateForChangedControls(model, state, forces)
            [model.parentModel, state] = model.parentModel.updateForChangedControls(state, forces);
        end

        function model = validateModel(model, varargin)
            model.parentModel = model.parentModel.validateModel(varargin{:});
            return
        end

        function [fn, index] = getVariableField(model, name)
            [fn, index] = model.parentModel.getVariableField(name);
        end

    end
    
end