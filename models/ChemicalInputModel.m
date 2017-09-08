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
        
        %%
        function model = ChemicalInputModel()
            model = model@ChemicalModel();
            model.inputNames = {};
        end
        
        %%
        function model = validateModel(model)
            model = validateModel@ChemicalModel(model);
            % setup unknownNames
            unknownNames = horzcat(model.CompNames, model.MasterCompNames, model.CombinationNames, model.SolidNames, 'poro'); % model.GasNames
            
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

            nPsi = sum(cellfun(@(x) ~isempty(x), regexpi( model.unknownNames, 'psi')));
%             assert(numel(model.unknownNames)-nPsi == model.nR + model.nMC + model.nLC + 1, 'well..., not as we are thinking...');

        end
        
        %%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            [pVars, logcomps, logmasterComps, comboComps, logGasComps, logSolidComps, logPoro]...
                = prepStateForEquations(model, state);

            [eqs, names, types] = equationsChemicalInit(logPoro, logcomps, logmasterComps, comboComps, ...
                                                       logGasComps, logSolidComps, state, model);

            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end

        %%
        function [logUnknowns, logComps, logMasterComps, combinationComps,...
                 logGasComps, logSolidComps, logPorosity] = prepStateForEquations(model, ...
                                                              state)
            CNames = model.logCompNames;
            MCNames = model.logMasterCompNames;
            LCNames = model.CombinationNames;
            GNames = model.logGasNames;
            SNames = model.logSolidNames;

            nC = numel(CNames);

            logUnknowns = model.logUnknownNames;
            logKnowns = model.logInputNames;

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
                mcInd = strcmpi(logUnknowns, GNames{i});
                if any(mcInd)
                    logGasComps{i} = logUnknownVal{mcInd};
                end
                mcInd = strcmpi(logKnowns, GNames{i});
                if any(mcInd)
                    logGasComps{i} = logKnownVal{mcInd};
                end
            end
            
            logSolidComps = cell(1,model.nS);
            for i = 1 : model.nS
                mcInd = strcmpi(logUnknowns, SNames{i});
                if any(mcInd)
                    logSolidComps{i} = logUnknownVal{mcInd};
                end
                mcInd = strcmpi(logKnowns, SNames{i});
                if any(mcInd)
                    logSolidComps{i} = logKnownVal{mcInd};
                end
            end
            
            pInd = strcmpi(logUnknowns, 'logporo');
            logPorosity = logUnknownVal{pInd}; 
            
        end
        
        %%
        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)

            [state, report] = updateState@ChemicalModel(model, state, problem, ...
                                                        dx, drivingForces);
            state = model.syncFromLog(state);
            
            LC = model.CombinationMatrix;
            CM = model.CompositionMatrix;

            for i = 1 : model.nLC
                combMaxMatrix = diag(1./LC(:, i));
                maxvals = combMaxMatrix*((state.combinationComponents)');
                maxvals = (min(maxvals))';
                state = model.capProperty(state, model.CombinationNames{i}, 1e-50, ...
                                          maxvals);
            end

%             for i = 1 : model.nC
%                 maxMatrix = diag(1./CM(:, i));
%                 maxvals = maxMatrix*((state.components)');
%                 maxvals = (min(maxvals))';
%                 state = model.capProperty(state, model.CompNames{i}, 1e-50, ...
%                                           maxvals);
%             end
           
            state = model.syncLog(state);
            
        end

        %%
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
