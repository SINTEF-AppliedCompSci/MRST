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
        function model = ChemicalInputModel(chemsys)
            model = model@ChemicalModel(chemsys);
            model.inputNames = chemsys.inputNames;
            unknownNames = horzcat(chemsys.speciesNames, ...
                                   chemsys.elementNames, ...
                                   chemsys.combinationNames, ...
                                   chemsys.solidNames, ...
                                   chemsys.gasNames, ...
                                   chemsys.surfaceActivityCoefficientNames);
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

            chemsys = model.chemicalSystem;

            CNames  = chemsys.logSpeciesNames;
            MCNames = chemsys.logElementNames;
            LCNames = chemsys.combinationNames;
            GNames  = chemsys.logGasNames;
            SNames  = chemsys.logSolidNames;
            SPNames = chemsys.logSurfaceActivityCoefficientNames;
            
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
