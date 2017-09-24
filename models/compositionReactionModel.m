classdef compositionReactionModel < ChemicalModel
    
    
    properties
        inputNames   % Names of variables that are used as inputs
        unknownNames % Names of variables that are used as inputs

    end
    
    methods
        %%
        function model = compositionReactionModel()
            model = model@ChemicalModel();
        end
        
        %%
        function model = validateModel(model)
            model = validateModel@ChemicalModel(model);
            unknownNames = horzcat(model.componentNames, model.masterComponentNames, model.combinationNames, model.solidNames, model.gasNames, 'fluidVolumeFraction');
            ind = ismember(unknownNames, model.inputNames);
            model.unknownNames = unknownNames(~ind);

        end
        
        %%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [pVars, logComponents, logMasterComponents, combinationComponents, logGasVolumeFractions, logSolidVolumeFractions, logFluidVolumeFraction] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsCompositionReactionGuess(model, state, logFluidVolumeFraction, logComponents, logMasterComponents, combinationComponents, logGasVolumeFractions, logSolidVolumeFractions);
            
            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end
        
        %%
        function [unknowns, logComponents, logMasterComponents, combinationComponents, logGasVolumeFractions, logSolidVolumeFractions, logFluidVolumeFraction] = prepStateForEquations(model, ...
                                                              state)
            
            CNames = model.logComponentNames;
            MCNames = model.logMasterComponentNames;
            LCNames = model.combinationNames;
            GNames = model.logGasNames;
            SNames = model.logSolidNames;
                        
            unknowns = model.unknownNames;
            knowns = model.inputNames;
            
            unknowns = addLogToNames(unknowns);
            knowns = addLogToNames(knowns);

            for i = 1 : model.nLC
                unknowns = regexprep(unknowns, ['log'  LCNames{i}],  LCNames{i});
                knowns = regexprep(knowns, ['log'  LCNames{i}],  LCNames{i});
            end
            
            unknownVal = cell(1,numel(unknowns));
            [unknownVal{:}] = model.getProps(state, unknowns{:});
            [unknownVal{:}] = initVariablesADI(unknownVal{:});
            
            knownVal = cell(1,numel(knowns));
            [knownVal{:}] = model.getProps(state, knowns{:});
            

            logComponents           = distributeVariable( CNames, knowns, unknowns, knownVal, unknownVal );
            logMasterComponents     = distributeVariable( MCNames, knowns, unknowns, knownVal, unknownVal );
            combinationComponents   = distributeVariable( LCNames, knowns, unknowns, knownVal, unknownVal );
            logGasVolumeFractions   = distributeVariable( GNames, knowns, unknowns, knownVal, unknownVal );
            logSolidVolumeFractions = distributeVariable( SNames, knowns, unknowns, knownVal, unknownVal );

            pInd = strcmpi(unknowns, 'logFluidVolumeFraction');
            logFluidVolumeFraction = unknownVal{pInd};   
            
        end
        
        %%
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

            solver = NonLinearSolver();
%             solver.maxIterations = 10;
%             solver.minIterations = 5;
            model.nonlinearTolerance = 1e-12;
            dt = 0; % dummy timestep
            drivingForces = []; % drivingForces;
            inputstate0 = inputstate;

            [state, failure, report] = solveMinistep(solver, model, inputstate, ...
                                                     inputstate0, dt, ...
                                                     drivingForces);

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
        
        if ~isempty(model.surfInfo)
            
            [state] = potentialGuess(model, state);
            state = model.syncFromLog(state);
        end
        
        report = [];
        end
    
    end

end
    