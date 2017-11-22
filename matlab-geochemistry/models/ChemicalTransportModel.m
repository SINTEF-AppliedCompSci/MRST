classdef ChemicalTransportModel < WaterModel

    properties
        % Chemical model for the chemistry
        chemicalModel
        % List of all the variable names for the chemistry
        chemical_fds;
        % List of all the variable names for the transport
        transport_fds;
        % Matrix to compute, for a given component, the amount that is
        % contained in the fluid (water)
        fluidMat;
        % Matrix to compute, for a given component, the amount that is
        % attached to the surface
        surfMat;
        plotIter % plot each iteration of the netwon solver
        plotFinal % plot final iteration of the newton solver


    end


    methods

        function model = ChemicalTransportModel(G, rock, fluid, chemicalLogModel, varargin)

            model = model@WaterModel(G, rock, fluid, varargin{:});
            model.chemicalModel = chemicalLogModel;
            model.chemical_fds = model.chemicalModel.getAllVarsNames();
            model.transport_fds = {'p', 'wellSol'};

            % Create a matrix of the components that are on surfaces
            chemModel = model.chemicalModel;
            nC        = chemModel.nC;
            nMC       = chemModel.nMC;
            CM        = chemModel.compositionMatrix;

            surfMaster  = logical(model.chemicalModel.surfMaster);
            surfComp    = sum(logical(CM(surfMaster, :)), 1);
            surfMult    = repmat(surfComp, nMC, 1);
            surfMatFlag = logical(CM.*surfMult);

            surfMat = zeros(nMC, nC);
            surfMat(surfMatFlag) = CM(surfMatFlag);
            model.surfMat = surfMat;
            
            fluidMat = zeros(nMC, nC);
            fluidMat(~surfMatFlag) = CM(~surfMatFlag);
            model.fluidMat = fluidMat;
            
            model.plotIter = false;
            model.plotFinal = false;


        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)

            [pVars, pressure, logComponents, logMasterComponents, combinations,...
                   logSaturationIndicies, logPartialPressures,...
                   logSurfaceAcitivityCoefficients] = prepStateForEquations(model, state);
               
            components = cellfun(@(x) exp(x), logComponents, 'UniformOutput',false);
            masterComponentss = cellfun(@(x) exp(x), logMasterComponents, 'UniformOutput',false);

            [chem_eqs, chem_names, chem_types] = equationsChemicalLog(model.chemicalModel, state, logComponents, logMasterComponents, combinations, ...
                                                       logPartialPressures, logSaturationIndicies,logSurfaceAcitivityCoefficients);


            [tr_eqs, tr_names, tr_types] = equationsTransportComponents(state0, ...
                                                              pressure, masterComponentss, ...
                                                              components,...
                                                              state, model, ...
                                                              dt, ...
                                                              drivingForces);
            eqs = horzcat(tr_eqs, chem_eqs );
            names = { tr_names{:},chem_names{:}};
            types = { tr_types{:},chem_types{:}};

            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end

        function [variableNames, pressure, logComponents, logMasterComponents, combinations,...
                   logSaturationIndicies, logPartialPressures,...
                   logSurfaceActivityCoefficients] = prepStateForEquations(model, state)
            
            chemModel = model.chemicalModel;

            logComponentNames       = chemModel.logSpeciesNames;
            logMasterComponentNames = chemModel.logElementNames;
            logSolidNames           = chemModel.logSolidNames;
            logGasNames             = chemModel.logGasNames;
            combinationNames        = chemModel.combinationNames;
            logSurfActNames         = chemModel.logSurfaceActivityCoefficientNames;
            
            
            variableNames = ['pressure', logComponentNames, logMasterComponentNames, logGasNames, logSolidNames, logSurfActNames];
            variableValues = cell(1, numel(variableNames));
            variableValues{1} = model.getProps(state, 'pressure');
            [variableValues{2:end}] = chemModel.getProps(state, variableNames{2:end});
            
            [variableValues{:}] = initVariablesADI(variableValues{:});
            
            logComponents        = cell(1, numel(logComponentNames));
            logMasterComponents  = cell(1, numel(logMasterComponentNames));
            logPartialPressures     = cell(1, numel(logGasNames));
            logSaturationIndicies   = cell(1, numel(logSolidNames));
            combinations   = cell(1, numel(combinationNames));
            logSurfaceActivityCoefficients = cell(1, numel(logSurfActNames));
            
            for i = 1 : numel(combinationNames)
                ind = strcmpi(combinationNames{i}, variableNames);
                combinations{i} = variableValues{ind};
            end

            for i = 1 : numel(logComponentNames)
                ind = strcmpi(logComponentNames{i}, variableNames);
                logComponents{i} = variableValues{ind};
            end

            for i = 1 : numel(logComponentNames)
                ind = strcmpi(logComponentNames{i}, variableNames);
                logComponents{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logMasterComponentNames)
                ind = strcmpi(logMasterComponentNames{i}, variableNames);
                logMasterComponents{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logGasNames)
                ind = strcmpi(logGasNames{i}, variableNames);
                logPartialPressures{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logSolidNames)
                ind = strcmpi(logSolidNames{i}, variableNames);
                logSaturationIndicies{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logSurfActNames)
                ind = strcmpi(logSurfActNames{i}, variableNames);
                logSurfaceActivityCoefficients{i} = variableValues{ind};
            end
            
            
            ind = strcmpi('pressure', variableNames);
            pressure = variableValues{ind};
      
               
               
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments

            chemModel = model.chemicalModel;

            vars = problem.primaryVariables;

            ind = false(size(vars));
            chemvars = {chemModel.logSpeciesNames{:}, chemModel.logElementNames{:}, chemModel.logGasNames{:}, chemModel.logSolidNames{:}, chemModel.logSurfaceActivityCoefficientNames{:}}; % the chemical primary variables, see getEquations
            [lia, loc] = ismember(chemvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;

            chem_problem                  = problem;
            chem_problem.primaryVariables = vars(ind);
            chem_dx                       = dx(ind);
            
%             state = chemModel.synclog(state);
            [state, chem_report] = chemModel.updateState(state, chem_problem, ...
                                                    chem_dx, drivingForces);


            ind = false(size(vars));
            fluidvars = {'pressure'}; % the chemical primary variables, see getEquations
            [lia, loc] = ismember(fluidvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;

            fluid_problem                  = problem;
            fluid_problem.primaryVariables = vars(ind);
            fluid_dx                       = dx(ind);

            [state, fluid_report] = updateState@WaterModel(model, state, ...
                                                           fluid_problem, ...
                                                           fluid_dx, ...
                                                           drivingForces);
            report = []; % no report for the moment.
            
            if model.plotIter
                h = findobj('tag', 'updateSpeciesfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updateSpeciesfig');
                    h = findobj('tag', 'updateSpeciesfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logSpecies);
                title('species - iteration');
                xlabel('cell index');
                ylabel('log_e species concentrations')
                legend(model.chemicalModel.speciesNames);
                drawnow;
                
                h = findobj('tag', 'updateMasterfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updateMasterfig');
                    h = findobj('tag', 'updateMasterfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logElements);
                title('elements - iteration');
                xlabel('cell index');
                ylabel('log_e element concentrations')
                legend(model.chemicalModel.elementNames);
                drawnow;
                
                
            end

        end

        function [state, report] = updateAfterConvergence(model, state0, state, ...
                                                          dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@WaterModel(model, state0, ...
                                                              state, dt, drivingForces);
           
            
            if model.plotFinal 
                h = findobj('tag', 'convergedfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'convergedfig');
                    h = findobj('tag', 'convergedfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logSpecies);
                title('components - converged');
                xlabel('cell index');
                ylabel('log_e species concentrations')
                legend(model.chemicalModel.speciesNames);

                h = findobj('tag', 'convergedmasterfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'convergedmasterfig');
                    h = findobj('tag', 'convergedmasterfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logElements);
                xlabel('cell index');
                ylabel('log_e element concentrations')
                title('master components - converged');
                legend(model.chemicalModel.elementNames);
                drawnow;
            end
            
        end


        function [eq, src] = addComponentContributions(model, cname, eq, ...
                                                       component, src, force)
            % Note: Here component denotes in fact the fluid part of the master component.
            if isempty(force)
                return
            end

            chemModel = model.chemicalModel;
            ind = strcmpi(cname, chemModel.elementNames);
            if chemModel.surfMaster(ind)
                return
            end

            c = model.getProp(force, cname);

            cells = src.sourceCells;
            qW = src.phaseMass{1}./model.fluid.rhoWS;

            isInj = qW > 0;

            qC = (isInj.*c + ~isInj.*component(cells)).*qW;

            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end

        function names = getComponentNames(model)
            names = model.chemicalModel.elementNames;
        end

        function [fn, index] = getVariableField(model, name)
            if ismember(name, model.chemical_fds)
                [fn, index] = model.chemicalModel.getVariableField(name);
            else
                [fn, index] = getVariableField@WaterModel(model, name);
            end
        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@WaterModel(model);
            % Sources for chemical components
            forces.chemSrc = [];
            % chemSrc is a struc with fields
            % 'cells' : injection cells
            % 'comps' : compositions (for each cell)
            % 'rates' : injection rates (for each cell)
        end

    end
end
