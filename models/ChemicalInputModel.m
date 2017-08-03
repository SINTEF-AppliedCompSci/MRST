classdef ChemicalInputModel < ChemicalModel
% model for initializing and solving the chemical system based on user chosen input

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

            solver = NonLinearSolver();
%             solver.maxIterations= 25;

            solver.LinearSolver.tolerance = 1e-14;
            model.nonlinearTolerance = 1e-14;
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
    