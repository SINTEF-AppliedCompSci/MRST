classdef compositionReactionModel < ChemicalModel
    
    
    properties
        inputNames   % Names of variables that are used as inputs
        unknownNames % Names of variables that are used as inputs
        logUnknownNames
        logInputNames % input names with log attached

    end
    
    methods

        function model = compositionReactionModel()
            model = model@ChemicalModel();
        end
        
        function model = validateModel(model)
            model = validateModel@ChemicalModel(model);
            % setup unknownNames
            unknownNames = horzcat(model.CompNames, model.MasterCompNames);
            ind = cellfun(@(name)(strcmpi(name, model.inputNames)), unknownNames, ...
                          'Uniformoutput', false);
            Pind = cellfun(@(x) ~isempty(x) , regexpi(unknownNames, 'psi'), 'Uniformoutput', false);
            
            Pind = cell2mat(Pind');
            ind = cell2mat(ind');
            ind = sum(ind, 2);
            ind = logical(ind) + Pind;
            model.unknownNames = unknownNames(~ind);
            model.logUnknownNames = cellfun(@(name) ['log', name], ...
                                            model.unknownNames, 'uniformoutput', ...
                                            false); 
            CNInd = cellfun(@(x) ~isempty(x) , regexpi(model.CompNames, 'psi'), 'Uniformoutput', false); 
            CNInd = cell2mat(CNInd');
            model.CompNames= model.CompNames(~CNInd);
            model.logCompNames = model.logCompNames(~CNInd);
            
            model.CompositionMatrix(:,CNInd) = [];
            model.ReactionMatrix(:,CNInd) = [];

            
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [logComps, logMasterComps] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsCompositionReactionGuess(logComps, logMasterComps, model);
            
            primaryVariables = model.logUnknownNames;
            problem = LinearizedProblem(eqs, types, names, primaryVariables, state, dt);

        end
        
        
        function [logComps, logMasterComps] = prepStateForEquations(model, ...
                                                              state)
            
            CNames = model.logCompNames;
            MCNames = model.logMasterCompNames;
            
            nC = numel(CNames);
            
            logUnknowns = model.logUnknownNames;
            logKnowns = model.logInputNames;

            logUnknownVal = cell(1,numel(logUnknowns));
            [logUnknownVal{:}] = model.getProps(state, logUnknowns{:});
            [logUnknownVal{:}] = initVariablesADI(logUnknownVal{:});
            
            logKnownVal = cell(1,numel(logKnowns));
            [logKnownVal{:}] = model.getProps(state, logKnowns{:});
            
            logComps = cell(1,nC);
            for i = 1 : nC
                cInd = strcmpi(logUnknowns, CNames{i});
                if any(cInd)
                    logComps{i} = logUnknownVal{cInd};
                end
                cInd = strcmpi(logKnowns, CNames{i});
                if any(cInd)
                    logComps{i} = logKnownVal{cInd};
                end
            end
            
            logMasterComps = cell(1,model.nMC);
            for i = 1 : model.nMC
                mcInd = strcmpi(logUnknowns, MCNames{i});
                if any(mcInd)
                    logMasterComps{i} = logUnknownVal{mcInd};
                end
                mcInd = strcmpi(logKnowns, MCNames{i});
                if any(mcInd)
                    logMasterComps{i} = logKnownVal{mcInd};
                end
            end
            
        end
        
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

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
            state = model.syncFromLog(state);
            
            nonLogVariables = regexprep(problem.primaryVariables, 'log', '');
            nC = numel(nonLogVariables);
            
            len = cellfun(@(x) length(x), nonLogVariables);
            [~,sortInd] = sort(len(:),1, 'ascend');
            pVar = nonLogVariables(sortInd);
            
            for i = 1 : nC
                
                p = pVar{i};
                compInd = strcmpi(p, model.CompNames);
                
                if any(strcmpi(p, model.MasterCompNames))
                    state = model.capProperty(state, p, eps, 2.5*mol/litre); 
                else
                    maxvals = model.maxMatrices{compInd}*((state.masterComponents)');
                    maxvals = (min(maxvals))';             
                    state = model.capProperty(state, p, eps, maxvals); 
                end
                
            end
            
            state = model.syncLog(state);

         end
        
        %
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
        % solve for the electrostatics
        
        if ~isempty(model.surfInfo)
            
            [state] = potentialGuess(model, state);
            state = model.syncFromLog(state);
        end
        
        report = [];
        end
    
    end

end
    