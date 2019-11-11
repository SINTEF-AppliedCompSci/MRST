classdef ChemicalInputModel < ChemicalModel
% Model for initializing and solving the chemical system based on user chosen
% input. See function initState of ChemicalModel.

%{
Copyright 2009-2016 SINTEF DIGITAL, Applied Mathematics and Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    properties
        inputNames   % Names of variables that are used as inputs
        unknownNames % Names of variables that are used as inputs
    end

    methods
        
        %%
        function model = ChemicalInputModel()
            model = model@ChemicalModel();
            model.inputNames = {};
        end
        
        %%
        function model = validateModel(model)
            model = validateModel@ChemicalModel(model);
            unknownNames = horzcat(model.speciesNames, model.elementNames, model.combinationNames, model.solidNames, model.gasNames, model.surfaceActivityCoefficientNames);
            ind = ismember(unknownNames, model.inputNames);
            model.unknownNames = unknownNames(~ind);

        end
        
        %%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [pVars, state] = prepStateForEquations(model, state);
            [eqs, names, types] = equationsChemicalLog(model, state);
            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end

        %%
        function [unknowns, state] = prepStateForEquations(model, state)
            
            CNames  = model.logSpeciesNames;
            MCNames = model.logElementNames;
            LCNames = model.combinationNames;
            GNames  = model.logGasNames;
            SNames  = model.logSolidNames;
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
            
            logSpecies                     = distributeVariable(CNames , knowns, unknowns, knownVal, unknownVal);
            logElements                    = distributeVariable(MCNames, knowns, unknowns, knownVal, unknownVal);
            combinationComponents          = distributeVariable(LCNames, knowns, unknowns, knownVal, unknownVal);
            logPartialPressures            = distributeVariable(GNames , knowns, unknowns, knownVal, unknownVal);
            logSaturationIndicies          = distributeVariable(SNames , knowns, unknowns, knownVal, unknownVal);
            logSurfaceActivityCoefficients = distributeVariable(SPNames, knowns, unknowns, knownVal, unknownVal);

            state = model.setProp(state, 'logSpecies', logSpecies);
            state = model.setProp(state, 'logElements', logElements);
            state = model.setProp(state, 'combinationComponents', combinationComponents);
            state = model.setProp(state, 'logPartialPressures', logPartialPressures);
            state = model.setProp(state, 'logSaturationIndicies', logSaturationIndicies);
            state = model.setProp(state, 'logSurfaceActivityCoefficients', logSurfaceActivityCoefficients);
            
        end

        
        %%
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

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



    end

end
