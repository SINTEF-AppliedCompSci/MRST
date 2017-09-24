classdef chargeBalanceModel < ChemicalInputModel
    
    
    properties
        CVC         % charge variation component
    end
    
    methods

        function model = chargeBalanceModel()
            model = model@ChemicalInputModel();
            model.CVC = [];
        end
        

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            if ~isfield(state, 'CVC')
            	nCells = size(state.components,1);
                state.CVC = zeros(nCells,1);
            end
            
            [unknowns, components, masterComponents, combinationComponents,...
                 gasVolumeFractions, solidVolumeFractions, fluidVolumeFraction, surfaceAcitivityCoefficients,CVC] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsChargeBalance(model, state, components, masterComponents, combinationComponents,...
                 gasVolumeFractions, solidVolumeFractions, fluidVolumeFraction, surfaceAcitivityCoefficients, CVC);
            
            problem = LinearizedProblem(eqs, types, names, unknowns, state, dt);

        end
        
        
        function [unknowns, components, masterComponents, combinationComponents,...
                 gasVolumeFractions, solidVolumeFractions, fluidVolumeFraction, surfaceAcitivityCoefficients,CVC] = prepStateForEquations(model, ...
                                                              state)
            
            CNames = model.logComponentNames;
            MCNames = model.logMasterComponentNames;
            LCNames = model.combinationNames;
            GNames = model.logGasNames;
            SNames = model.logSolidNames;
            SPNames = model.logSurfaceActivityCoefficientNames;
            
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
            

            components           = distributeVariable( CNames, knowns, unknowns, knownVal, unknownVal );
            masterComponents     = distributeVariable( MCNames, knowns, unknowns, knownVal, unknownVal );
            combinationComponents   = distributeVariable( LCNames, knowns, unknowns, knownVal, unknownVal );
            gasVolumeFractions   = distributeVariable( GNames, knowns, unknowns, knownVal, unknownVal );
            solidVolumeFractions = distributeVariable( SNames, knowns, unknowns, knownVal, unknownVal );
            surfaceAcitivityCoefficients = distributeVariable( SPNames, knowns, unknowns, knownVal, unknownVal );

            pInd = strcmpi(unknowns, 'logFluidVolumeFraction');
            fluidVolumeFraction = unknownVal{pInd}; 

            cvcInd = strcmpi(unknowns, 'CVC');
            CVC = unknownVal{cvcInd}; 
            
        end
        
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            % grab the names of unknowns                                              
            unknownNames = horzcat(model.componentNames, model.masterComponentNames,model.combinationNames, 'CVC', 'fluidVolumeFraction', model.solidNames, model.gasNames, model.surfaceActivityCoefficientNames);
            ind = ismember(unknownNames, model.inputNames);
            model.unknownNames = unknownNames(~ind);
            
            solver = NonLinearSolver();   
            solver.LinearSolver.tolerance = 1e-12;
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
            state0 = state;
            [state, report] = updateState@PhysicalModel(model, state, problem, ...
                                                        dx, drivingForces);
        
            state = model.syncFromLog(state);                                       
            state = updateChemicalModel(model, problem, state, state0 );
            state = model.syncLog(state);
             
                                                   
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces) %#ok
        % solve for the electrostatics
        
            CVC = model.getProp(state, 'CVC');
            MCval = model.getProp(state, model.CVC);
            
            warningText = ['Charge balance could not be acheived with given constraint on ' model.CVC ' increase the value if reasonable.'];
            
            if ~all(abs(CVC) <= 2*MCval),
                warning(warningText);
            end
            
            state = model.setProp(state, model.CVC, MCval + CVC);
            assert(all((MCval + CVC) > 0), warningText);
           
            
            state = model.syncLog(state);
        
        
            report = [];
        end
    
    end

end
    