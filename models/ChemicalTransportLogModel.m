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


        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            [p, wellSol] = model.getProps(state, 'pressure', 'wellSol');


            chemModel = model.chemicalModel;
            % chemical equations

            compNames = chemModel.CompNames;
            masterCompNames = chemModel.MasterCompNames;
            logCompNames = chemModel.logCompNames;
            logMasterCompNames = chemModel.logMasterCompNames;

            comps = cell(1, numel(compNames));
            mastercomps = cell(1, numel(masterCompNames));

            logComps = cell(1, numel(logCompNames));
            logMasterComps = cell(1, numel(logMasterCompNames));

            [logComps{:}] = model.getProps(state, logCompNames{:});
            [logMasterComps{:}] = model.getProps(state, logMasterCompNames{:});

            [p, logComps{:}, logMasterComps{:}] = initVariablesADI(p, logComps{:}, ...
                                                             logMasterComps{:});

            comps = cellfun(@(x) exp(x), logComps, 'UniformOutput',false);
            masterComps = cellfun(@(x) exp(x), logMasterComps, 'UniformOutput',false);

            [chem_eqs, chem_names, chem_types] = equationsChemicalLog(logComps, ...
                                                              logMasterComps, ...
                                                              chemModel);


            [tr_eqs, tr_names, tr_types] = equationsTransportComponents(state0, ...
                                                              p, masterComps, ...
                                                              comps,...
                                                              state, model, ...
                                                              dt, ...
                                                              drivingForces);

            primaryVars = {'pressure', chemModel.logCompNames{:}, chemModel.logMasterCompNames{:}};
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
            chemvars = {chemModel.logCompNames{:}, chemModel.logMasterCompNames{:}}; % the chemical primary variables, see getEquations
            [lia, loc] = ismember(chemvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;

            chem_problem                  = problem;
            chem_problem.primaryVariables = vars(ind);
            chem_dx                       = dx(ind);

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

        function [state, report] = updateAfterConvergence(model, state0, state, ...
                                                          dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@WaterModel(model, state0, ...
                                                              state, dt, drivingForces);

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
