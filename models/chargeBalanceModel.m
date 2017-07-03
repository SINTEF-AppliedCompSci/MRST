classdef chargeBalanceModel < ChemicalInputModel
    
    
    properties
        CVC         % charge variation component
    end
    
    methods

        function model = chargeBalanceModel()
            model = model@ChemicalInputModel();
            model.CVC = [];
        end
        
%         function model = validateModel(model)
% %             model = validateModel@ChemicalInputModel(model);
%             % setup unknownNames
%             unknownNames = horzcat(model.CompNames, model.MasterCompNames, 'CVC');
%             ind = cellfun(@(name)(strcmpi(name, model.inputNames)), unknownNames, ...
%                           'Uniformoutput', false);
% 
%             model.unknownNames = unknownNames(~ind);
%             
%             
%         end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            if size(state.components,2) ~= numel(model.unknownNames)
            	nCells = size(state.components,1);
                state.components(:,end+1) = zeros(nCells,1);
            end
            
            [comps, masterComps] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsChargeBalance(model, comps, masterComps);
            
            primaryVariables =model.unknownNames;
            problem = LinearizedProblem(eqs, types, names, primaryVariables, state, dt);

        end
        
        
        function [comps, masterComps] = prepStateForEquations(model, ...
                                                              state)
            
            CNames = horzcat(model.CompNames);
            MCNames = model.MasterCompNames;
            unknowns = horzcat(model.unknownNames);
            knowns = model.inputNames;
            
            nC = numel(CNames);
            
            
            unknownVal = cell(1,numel(unknowns));
            [unknownVal{:}] = model.getProps(state, unknowns{:});
            [unknownVal{:}] = initVariablesADI(unknownVal{:});
            
            knownVal = cell(1,numel(knowns));
            [knownVal{:}] = model.getProps(state, knowns{:});
            
            comps = cell(1,nC);
            for i = 1 : nC
                cInd = strcmpi(unknowns, CNames{i});
                if any(cInd)
                    comps{i} = unknownVal{cInd};
                end
                cInd = strcmpi(knowns, CNames{i});
                if any(cInd)
                    comps{i} = knownVal{cInd};
                end
            end
            
            masterComps = cell(1,model.nMC);
            for i = 1 : model.nMC
                mcInd = strcmpi(unknowns, MCNames{i});
                if any(mcInd)
                    masterComps{i} = unknownVal{mcInd};
                end
                mcInd = strcmpi(knowns, MCNames{i});
                if any(mcInd)
                    masterComps{i} = knownVal{mcInd};
                end
            end

        end
        
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            % grab the names of unknowns                                              
            unknownNames = horzcat(model.CompNames, model.MasterCompNames, 'CVC');
            ind = cellfun(@(name) strcmpi(name, model.inputNames), unknownNames, ...
                          'Uniformoutput', false);
                      
            ind = cell2mat(ind');
            ind = sum(ind, 2);

            model.unknownNames = unknownNames(~ind);
            model.CompNames = horzcat(model.CompNames, 'CVC');
            
            solver = NonLinearSolver();
            model.nonlinearTolerance = 1e-12;
            dt = 0; % dummy timestep
            drivingForces = []; % drivingForces;
            inputstate0 = inputstate;

            [state, failure, report] = solveMinistep(solver, model, inputstate, ...
                                                     inputstate0, dt, ...
                                                     drivingForces);

        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments
            [state, report] = updateState@PhysicalModel(model, state, problem, ...
                                                        dx, drivingForces);
                                                    
            for i = 1 : numel(problem.primaryVariables)
                p = problem.primaryVariables{i};
                if ~strcmpi(p, 'CVC')
                    state = model.capProperty(state, p, eps);
                end
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
        % solve for the electrostatics
        
            CVC = model.getProp(state, 'CVC');
            MCval = model.getProp(state, model.CVC);
            
            if ~all(abs(CVC) <= 2*MCval),
                warning(['Charge balance could not be acheived with given constraint on ' model.CVC ' increase the value if reasonable.']);
            end
            
            state = model.setProp(state, model.CVC, MCval + CVC);
            assert(all((MCval + CVC) > 0), ['Charge balance could not be acheived with given constraint on ' model.CVC ' increase the value if reasonable. Or decrease the concentration of other charged species in the system.']);

            state.components(:,end) = [];
           
            
            state = model.syncLog(state);
        
        
            report = [];
        end
    
    end

end
    