classdef ChemicalTransportLogModel < WaterModel

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
        plotIter %to plot or not to plot

    end


    methods

        function model = ChemicalTransportLogModel(G, rock, fluid, chemicalLogModel, varargin)

            model = model@WaterModel(G, rock, fluid, varargin{:});
            model.chemicalModel = chemicalLogModel;
            model.chemical_fds = model.chemicalModel.getAllVarsNames();
            model.transport_fds = {'p', 'wellSol'};

            % Create a matrix of the components that are on surfaces
            chemModel = model.chemicalModel;
            nC        = chemModel.nC;
            nMC       = chemModel.nMC;
            CM        = chemModel.CompositionMatrix;

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


        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)


            chemModel = model.chemicalModel;
            % chemical equations

            logCompNames = chemModel.logCompNames;
            logMasterCompNames = chemModel.logMasterCompNames;
            logSolidCompNames = chemModel.logSolidNames;
            logGasCompNames = chemModel.logGasNames;
            
            variableNames = ['pressure', logCompNames, logMasterCompNames, logGasCompNames, logSolidCompNames, 'logPoro'];
            variableValues = cell(1, numel(variableNames));
            variableValues{1} = model.getProps(state, 'pressure');
            [variableValues{2:end}] = chemModel.getProps(state, variableNames{2:end});
            
            [variableValues{:}] = initVariablesADI(variableValues{:});
            
            logComps        = cell(1, numel(logCompNames));
            logMasterComps  = cell(1, numel(logMasterCompNames));
            logGasComps     = cell(1, numel(logGasCompNames));
            logSolidComps   = cell(1, numel(logSolidCompNames));
            
            for i = 1 : numel(logCompNames)
                ind = strcmpi(logCompNames{i}, variableNames);
                logComps{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logMasterCompNames)
                ind = strcmpi(logMasterCompNames{i}, variableNames);
                logMasterComps{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logGasCompNames)
                ind = strcmpi(logGasCompNames{i}, variableNames);
                logGasComps{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logSolidCompNames)
                ind = strcmpi(logSolidCompNames{i}, variableNames);
                logSolidComps{i} = variableValues{ind};
            end
            
            ind = strcmpi('logPoro', variableNames);
            logPoro = variableValues{ind};
            
            ind = strcmpi('pressure', variableNames);
            p = variableValues{ind};

            comps = cellfun(@(x) exp(x), logComps, 'UniformOutput',false);
            masterComps = cellfun(@(x) exp(x), logMasterComps, 'UniformOutput',false);

            [chem_eqs, chem_names, chem_types] = equationsChemicalLog(logPoro, logComps,...
                                                    logMasterComps, logGasComps,...
                                                    logSolidComps, state, chemModel);


            [tr_eqs, tr_names, tr_types] = equationsTransportComponents(state0, ...
                                                              p, masterComps, ...
                                                              comps,logPoro,...
                                                              state, model, ...
                                                              dt, ...
                                                              drivingForces);

            primaryVars = {'pressure', logCompNames{:}, logMasterCompNames{:}, logSolidCompNames{:}, logGasCompNames{:}, 'logPoro'};
            eqs = horzcat(tr_eqs, chem_eqs );
            names = { tr_names{:},chem_names{:}};
            types = { tr_types{:},chem_types{:}};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end


        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments

            chemModel = model.chemicalModel;

            vars = problem.primaryVariables;

            ind = false(size(vars));
            chemvars = {chemModel.logCompNames{:}, chemModel.logMasterCompNames{:}, chemModel.logGasNames{:}, chemModel.logSolidNames{:}, 'logPoro'}; % the chemical primary variables, see getEquations
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
                h = findobj('tag', 'updatefig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updatefig');
                    h = findobj('tag', 'updatefig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(log10(state.components*litre/mol));
                title('components');
                legend(model.chemicalModel.CompNames);

                drawnow;
            end

        end

        function [state, report] = updateAfterConvergence(model, state0, state, ...
                                                          dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@WaterModel(model, state0, ...
                                                              state, dt, drivingForces);
           
            rPoro = model.rock.poro;      
            stepPoro = [state.poro state.solidComponents];
            for i = 1 : size(stepPoro,2);
                stepPoro(:,i) = stepPoro(:,i).*rPoro;
            end
            
            h = findobj('tag', 'convergedPorofig');
            if isempty(h)
                figure
                set(gcf, 'tag', 'convergedPorofig');
                h = findobj('tag', 'convergedPorofig');
            end
            set(0, 'currentfigure', h)
            clf
            plot(stepPoro);
            title('porosities - converged');
            legend(['porosity' model.chemicalModel.SolidNames]);
            
            
            h = findobj('tag', 'convergedfig');
            if isempty(h)
                figure
                set(gcf, 'tag', 'convergedfig');
                h = findobj('tag', 'convergedfig');
            end
            set(0, 'currentfigure', h)
            clf
            plot(log10(state.components*litre/mol));
            title('components - converged');
            legend(model.chemicalModel.CompNames);

            h = findobj('tag', 'convergedmasterfig');
            if isempty(h)
                figure
                set(gcf, 'tag', 'convergedmasterfig');
                h = findobj('tag', 'convergedmasterfig');
            end
            set(0, 'currentfigure', h)
            clf
            plot(log10(state.masterComponents*litre/mol));
            title('master components - converged');
            legend(model.chemicalModel.MasterCompNames);
            drawnow;

        end


        function [eq, src] = addComponentContributions(model, cname, eq, ...
                                                       component, src, force)
            % Note: Here component denotes in fact the fluid part of the master component.
            if isempty(force)
                return
            end

            chemModel = model.chemicalModel;
            ind = strcmpi(cname, chemModel.MasterCompNames);
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
            names = model.chemicalModel.MasterCompNames;
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
