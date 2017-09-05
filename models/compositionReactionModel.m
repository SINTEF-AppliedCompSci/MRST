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
            unknownNames = horzcat(model.CompNames, model.MasterCompNames, model.CombinationNames, model.GasNames, model.SolidNames, 'poro');
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

            [pVars, logComps, logMasterComps, comboComps, logGasComps, logSolidComps, logPorosity] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsCompositionReactionGuess(state, logPorosity, logComps, logMasterComps, comboComps, logGasComps, logSolidComps, model);
            
            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end
        
        
        function [logUnknowns, logComps, logMasterComps, combinationComps, logGasComps, logSolidComps, logPorosity] = prepStateForEquations(model, ...
                                                              state)
            
            CNames = model.logCompNames;
            MCNames = model.logMasterCompNames;
            LCNames = model.CombinationNames;
            GNames = model.logGasNames;
            SNames = model.logSolidNames;
            
            nC = numel(CNames);
            
            logUnknowns = model.logUnknownNames;
            logKnowns = model.logInputNames;

            % actually, we want the non log form of the linear combination
            % variables
            for i = 1 : model.nLC
                logUnknowns = regexprep(logUnknowns, ['log'  LCNames{i}],  LCNames{i});
                logKnowns = regexprep(logKnowns, ['log'  LCNames{i}],  LCNames{i});
            end
            
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
            
            combinationComps = cell(1,model.nLC);
            for i = 1 : model.nLC
                mcInd = strcmpi(logUnknowns, LCNames{i});
                if any(mcInd)
                    combinationComps{i} = logUnknownVal{mcInd};
                end
                mcInd = strcmpi(logKnowns, LCNames{i});
                if any(mcInd)
                    combinationComps{i} = logKnownVal{mcInd};
                end
            end
            
            logGasComps = cell(1,model.nG);
            for i = 1 : model.nG
                gInd = strcmpi(logUnknowns, GNames{i});
                if any(gInd)
                    logGasComps{i} = logUnknownVal{gInd};
                end
                gInd = strcmpi(logKnowns, GNames{i});
                if any(gInd)
                    logGasComps{i} = logKnownVal{gInd};
                end
            end
           
            logSolidComps = cell(1,model.nS);
            for i = 1 : model.nS
                sInd = strcmpi(logUnknowns, SNames{i});
                if any(sInd)
                    logSolidComps{i} = logUnknownVal{sInd};
                end
                sInd = strcmpi(logKnowns, SNames{i});
                if any(sInd)
                    logSolidComps{i} = logKnownVal{sInd};
                end
            end

            
            pInd = strcmpi(logUnknowns, 'logporo');
            logPorosity = logUnknownVal{pInd};   
            
        end
        
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

            solver = NonLinearSolver();
            solver.maxIterations = 10;
            solver.minIterations = 5;
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
                elseif ismember(p, model.CombinationNames)
                    state = model.capProperty(state, p, -2.5*mol/litre, 2.5*mol/litre);
                elseif ismember(p, [model.GasNames, model.SolidNames, 'poro'])
                    state = model.capProperty(state, p, 1e-50, 1);
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
    