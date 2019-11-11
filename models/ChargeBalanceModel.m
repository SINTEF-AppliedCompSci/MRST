classdef ChargeBalanceModel < ChemicalInputModel
    
    
    properties
    end
    
    methods

        function model = ChargeBalanceModel()
            model = model@ChemicalInputModel();
        end
        

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            [unknowns, state] = prepStateForEquations(model, state);
            [eqs, names, types] = equationsChargeBalance(model, state);
            problem = LinearizedProblem(eqs, types, names, unknowns, state, dt);
            
        end

        function [unknowns, state] = prepStateForEquations(model, state)
                
            CNames  = model.logSpeciesNames;
            MCNames = model.logElementNames;
            LCNames = model.combinationNames;
            GNames  = model.logGasNames;
            SNames  = model.logSolidNames;
            SPNames = model.logSurfaceActivityCoefficientNames;
            
            unknowns = model.unknownNames;
            knowns = model.inputNames(~strcmpi(model.inputNames,model.CVC));
            
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
            

            logSpecies                      = distributeVariable(CNames , knowns, unknowns, knownVal, unknownVal);
            logElements                     = distributeVariable(MCNames, knowns, unknowns, knownVal, unknownVal);
            combinationComponents        = distributeVariable(LCNames, knowns, unknowns, knownVal, unknownVal);
            logPartialPressures             = distributeVariable(GNames , knowns, unknowns, knownVal, unknownVal);
            logSaturationIndicies          = distributeVariable(SNames , knowns, unknowns, knownVal, unknownVal);
            logSurfaceActivityCoefficients = distributeVariable(SPNames, knowns, unknowns, knownVal, unknownVal);

            state = model.setProp(state, 'logSpecies', logSpecies);
            state = model.setProp(state, 'logElements', logElements);
            state = model.setProp(state, 'combinationComponents', combinationComponents);
            state = model.setProp(state, 'logPartialPressures', logPartialPressures);
            state = model.setProp(state, 'logSaturationIndicies', logSaturationIndicies);
            state = model.setProp(state, 'logSurfaceActivityCoefficients', logSurfaceActivityCoefficients);
            
        end
        
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            % grab the names of unknowns                                              
            unknownNames = horzcat(model.speciesNames, ...
                                   model.elementNames, ...
                                   model.combinationNames, ...
                                   model.solidNames, ...
                                   model.gasNames, ...
                                   model.surfaceActivityCoefficientNames);
            
            ind = ismember(unknownNames, model.inputNames(~strcmpi(model.inputNames, model.CVC)));
            model.unknownNames = unknownNames(~ind);
            
            solver = NonLinearSolver(); 
            
            solver.maxIterations= model.nonLinearMaxIterations;
            solver.minIterations= model.nonLinearMinIterations;
            solver.LinearSolver.tolerance = model.linearTolerance;
            solver.LinearSolver.maxIterations = model.linearMaxIterations;
            
            dt = 0; % dummy timestep
            drivingForces = []; % drivingForces;
            inputstate0 = inputstate;
            
            [state, failure, report] = solveMinistep(solver, model, inputstate, ...
                                                     inputstate0, dt, ...
                                                     drivingForces);

        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments
            state0 = state;
            [state, report] = updateState@PhysicalModel(model, state, problem, ...
                                                        dx, drivingForces);
        
            state = model.syncFromLog(state);                                       
            state = updateChemicalModel(model, problem, state, state0 );
            state = model.syncLog(state);
                                                   
        end
        
    
    end

end
    