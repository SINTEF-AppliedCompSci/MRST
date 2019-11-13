classdef ChargeBalanceModel < ChemicalInputModel
    
    properties
        CVC % Charge variation component
    end
    
    methods

        function model = ChargeBalanceModel(chemsys)
            model = model@ChemicalInputModel(chemsys);
        end
        

        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
            [unknowns, state] = prepStateForEquations(model, state);
            [chem_eqs, chem_names, chem_types]       = equationsChemicalLog(model, state);
            [charge_eqs, charge_names, charge_types] = equationsChargeBalance(model, state);
            % we concatenate the equations
            eqs   = horzcat(chem_eqs, charge_eqs);
            names = {chem_names{:}, charge_names{:}};
            types = {chem_types{:}, charge_types{:}};
            
            problem = LinearizedProblem(eqs, types, names, unknowns, state, dt);
            
        end

        function [unknowns, state] = prepStateForEquations(model, state)
                
            chemsys = model.chemicalSystem;
            
            unknowns = model.unknownNames;
            knowns = model.inputNames(~strcmpi(model.inputNames, model.CVC));
            
            unknowns = addLogToNames(unknowns);
            knowns = addLogToNames(knowns);
            
            for i = 1 : chemsys.nLC
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
        
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            chemsys = model.chemicalSystem;
            % grab the names of unknowns                                              
            unknownNames = horzcat(chemsys.speciesNames, ...
                                   chemsys.elementNames, ...
                                   chemsys.combinationNames, ...
                                   chemsys.solidNames, ...
                                   chemsys.gasNames, ...
                                   chemsys.surfaceActivityCoefficientNames);
            
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
    