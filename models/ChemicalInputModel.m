classdef ChemicalInputModel < ChemicalModel
    
    
    properties
        inputNames   % Names of variables that are used as inputs
        unknownNames % Names of variables that are used as inputs
        logUnknownNames % names with log attached
        logInputNames % input names with log attached
    end
    
    methods

        function model = ChemicalInputModel()
            model = model@ChemicalModel();
            model.inputNames = {};
        end
        
        function model = validateModel(model)
            model = validateModel@ChemicalModel(model);
            % setup unknownNames
            unknownNames = horzcat(model.CompNames, model.MasterCompNames);
            ind = cellfun(@(name)(strcmpi(name, model.inputNames)), unknownNames, ...
                          'Uniformoutput', false);
            ind = cell2mat(ind');
            ind = sum(ind, 2);
            model.unknownNames = unknownNames(~ind);
            
            model.logUnknownNames = cellfun(@(name) ['log', name], ...
                                            model.unknownNames, 'uniformoutput', ...
                                            false);
                                        
            model.logInputNames = cellfun(@(name) ['log', name], ...
                                            model.inputNames, 'uniformoutput', ...
                                            false);                         
                                           
            assert(numel(model.unknownNames) == model.nC, 'well..., not as we are thinking...');
            
        end

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [logcomps, logmasterComps] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsChemicalLog(logcomps, logmasterComps, ...
                                                        model);
            
            primaryVariables = model.logUnknownNames;
            problem = LinearizedProblem(eqs, types, names, primaryVariables, state, dt);

        end
        
        function [logComps, logMasterComps] = prepStateForEquations(model, ...
                                                              state)
            
            logCompNames       = model.logCompNames;
            logMasterCompNames = model.logMasterCompNames;
            logUnknownNames    = model.logUnknownNames;
             
            logComps = cell(model.nC, 1);
            [logComps{:}] = model.getProps(state, model.logCompNames{:});

            logMasterComps = cell(model.nMC, 1);
            [logMasterComps{:}] = model.getProps(state, model.logMasterCompNames{:});

            logvals = cell(model.nC, 1);
            [logvals{:}] = model.getProps(state, model.logUnknownNames{:});
            [logvals{:}] = initVariablesADI(logvals{:});
            
            for i = 1 : model.nC
                indComp = strcmp(logUnknownNames{i}, logCompNames);
                indMComp = strcmp(logUnknownNames{i}, logMasterCompNames);
                if any(indComp)
                    logComps{indComp} = logvals{i};
                end
                if any(indMComp)
                    logMasterComps{indMComp} = logvals{i};
                end
            end
            

        end
        
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

            solver = NonLinearSolver('maxIterations', 100);
            dt = 0; % dummy timestep
            drivingForces = []; % drivingForces;
            inputstate0 = inputstate;

            [state, failure, report] = solveMinistep(solver, model, inputstate, ...
                                                     inputstate0, dt, ...
                                                     drivingForces);
            tmp = cell(1,model.nMC);
            [tmp{:}] = model.getProps(state, model.MasterCompNames{:});
            state.masterComponents  = horzcat(tmp{:});

            tmp = cell(1,model.nC);
            [tmp{:}] = model.getProps(state, model.CompNames{:});
            state.components        =  horzcat(tmp{:});
        end
        
       
    
    end

end
    