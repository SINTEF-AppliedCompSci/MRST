classdef CompositionReactionModel < ChemicalModel
    
    
    properties
        inputNames   % Names of variables that are used as inputs
        unknownNames % Names of variables that are used as inputs

    end
    
    methods
        %%
        function model = CompositionReactionModel(chemsys)
            model = model@ChemicalModel(chemsys);
            inputNames = chemsys.inputNames;
            
            unknownNames = horzcat(chemsys.speciesNames, ...
                                   chemsys.elementNames, ...
                                   chemsys.combinationNames, ...
                                   chemsys.solidNames, ...
                                   chemsys.gasNames);
            ind = ismember(unknownNames, inputNames);
            unknownNames = unknownNames(~ind);
            
            model.inputNames   = inputNames;
            model.unknownNames = unknownNames;
            
        end

        
        %%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [pVars, state] = prepStateForEquations(model, state);
            [eqs, names, types] = equationsCompositionReactionGuess(model, state);
            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end

        function [unknowns, state] = prepStateForEquations(model, state)
            
            chemsys = model.chemicalSystem;
                        
            unknowns = model.unknownNames;
            knowns   = model.inputNames;
            
            unknowns = addLogToNames(unknowns);
            knowns = addLogToNames(knowns);

            for i = 1 : chemsys.nLC
                % we do not use the logarithmic variables for the linear
                % combinations (they do not have necessarily a positive sign)
                LCNames = chemsys.combinationNames;
                unknowns = regexprep(unknowns, ['log'  LCNames{i}],  LCNames{i});
                knowns = regexprep(knowns, ['log'  LCNames{i}],  LCNames{i});
            end
            
            unknownVal = cell(1,numel(unknowns));
            [unknownVal{:}] = model.getProps(state, unknowns{:});
            [unknownVal{:}] = initVariablesADI(unknownVal{:});
            
            knownVal = cell(1,numel(knowns));
            [knownVal{:}] = model.getProps(state, knowns{:});
            
            names = {unknowns{:}, knowns{:}};
            vals  = {unknownVal{:}, knownVal{:}};

            state = model.setProps(state, names, vals);
            
        end
        
        %%
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            chemsys = model.chemicalSystem;
            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

            solver = NonLinearSolver();
            dt = 0; % dummy timestep
            drivingForces = []; % drivingForces;
            inputstate0 = inputstate;
%             solver.maxIterations = 10;

            [state, failure, report] = solveMinistep(solver, model, inputstate, ...
                                                     inputstate0, dt, ...
                                                     drivingForces);
            if ~isempty(chemsys.surfInfo)

                [state] = potentialGuess(model, state);
                state = model.syncFromLog(state);
            end
        end
        
        %%
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
            % Update state based on Newton increments
        	[state, report] = updateState@ChemicalModel(model, state, problem, ...
                                                        dx, drivingForces);
             state = model.syncLog(state);

         end
        
        %%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
        % solve for the electrostatics
        

        
        report = [];
        end
    
    end

end
    